import scipy.stats
import scipy.special
import numpy as np

import dataproducts

class RatingError(Exception):
    pass

class RatingDepreciatedError(RatingError):
    """ Error to raise when a rating version has been 
        superceded by a more recent version.""" 
    pass

class RatingWarning(Warning):
    pass


class RaterLoadWarning(RatingWarning):
    pass


class RatingInstanceIDCache(object):
    """A cache object to keep track of rating instance ID numbers.
        This object will query the database for ID numbers if
        they are not already cached. If not ID number is found
        in the database the cache obeject will add an entry
        to the DB.
    """
    def __init__(self, dbname='default'):
        """Constructor for RatingIDCache objects.

            Input:
                dbname: The database to connect to.
                    (Default: common-copy).

            Output:
                cache: The RatingIDCache object.
        """
        # Put import here so testing on systems 
        # without all prereqs doesn't fail.
        global database
        import database 
        self.dbname = dbname
        self.sp_idcache = {} # cache for SP rating instances
        self.pdm_idcache = {} # cache for PDM rating instances

        self.get_id = self.get_pdm_id # for backwards compatability

    def get_pdm_id(self, name, version, description):
        """Return the pdm_rating_instance_id of the rating
            matching the arguments provided. If there is
            no match cached fetch from the DB. If there is
            no match in the DB insert it. Any newly added/fetched
            IDs are added to the cache.
            
            Inputs:
                name: The rating's name.
                version: The rating's version
                description: The rating's description

            Output:
                id: The pdm_rating_type_id from the DB.
        """
        if (name, version) in self.pdm_idcache:
            id = self.pdm_idcache[(name,version)]
        else:
            # Check database for a match
            id = self._get_id_from_db(name, version, description, pdm_or_sp='pdm')
            self.pdm_idcache[(name, version)] = id
        return id

    def get_sp_id(self, name, version, description):
        """Return the sp_rating_instance_id of the rating
            matching the arguments provided. If there is
            no match cached fetch from the DB. If there is
            no match in the DB insert it. Any newly added/fetched
            IDs are added to the cache.
            
            Inputs:
                name: The rating's name.
                version: The rating's version
                description: The rating's description

            Output:
                id: The sp_rating_type_id from the DB.
        """
        if (name, version) in self.sp_idcache:
            id = self.sp_idcache[(name,version)]
        else:
            # Check database for a match
            id = self._get_id_from_db(name, version, description, pdm_or_sp='sp')
            self.sp_idcache[(name, version)] = id
        return id

    def _get_id_from_db(self, name, version, description, pdm_or_sp='pdm'):
        """Get the rating instance from the DB.

            Inputs:
                name: The rating's name.
                version: The rating's version
                description: The rating's description
                pdm_or_sp: Whether the rating is for PDM or SP candidates

            Output:
                id: The pdm_rating_type_id or sp_rating_type_id from the DB.
        """
        db = database.Database(self.dbname)
        query = "SELECT ri.%s_rating_instance_id, " \
                        "rt.name, " \
                        "ri.version " \
                   "FROM %s_rating_instance AS ri WITH(NOLOCK) " \
                   "LEFT JOIN %s_rating_type AS rt WITH(NOLOCK) " \
                        "ON rt.%s_rating_type_id=ri.%s_rating_type_id " \
                   "WHERE rt.name=? AND ri.version>=? " \
                   "ORDER BY ri.version DESC" % tuple([pdm_or_sp]*5)
        db.execute(query,(name, version))
        row = db.cursor.fetchone()
        db.close()
        if (row is not None) and (row[2] == version):
            id = row[0]
        elif (row is not None) and (row[2] > version):
            # We're considering a rating version that has
            # been superceded! Error...
            raise RatingDepreciatedError("This rating (%s (v%d)) is no longer " \
                              "current. It has been replaced by " \
                              "%s (v%d)" % \
                              (name, version, row[1], row[2]))
        else: # no row or row[2] < version
            id = self._add_id_to_db(name, version, description,pdm_or_sp=pdm_or_sp)
        return id

    def _add_id_to_db(self, name, version, description, pdm_or_sp="pdm"):
        """Add the rating instance (and type, if neccessary)
            to the database, and return the corresponding 
            pdm_rating_instance_id.

            Inputs:
                name: The rating's name.
                version: The rating's version
                description: The rating's description
                pdm_or_sp: Whether the rating is for PDM or SP candidates

            Output:
                id: The pdm_rating_type_id or sp_rating_type_id from the DB.
        """
        db = database.Database(self.dbname)
        # Check for rating type
        query = "SELECT %s_rating_type_id FROM %s_rating_type " \
                    "WHERE name=?" % tuple(2*[pdm_or_sp])
        db.execute(query, name)
        row = db.cursor.fetchone()
        
        # Add rating type if necessary
        if row is not None:
            type_id = row[0]
        else:
            query = "INSERT INTO %s_rating_type " \
                     "(name, description) " \
                     "VALUES (?, ?)" % pdm_or_sp
            db.execute(query, (name, description))
            query = "SELECT %s_rating_type_id " \
                        "FROM %s_rating_type " \
                        "WHERE name=?" % tuple(2*[pdm_or_sp])
            db.execute(query, name)
            type_id = db.cursor.fetchone()[0]

        # Check for rating instance
        query = "SELECT %s_rating_instance_id FROM %s_rating_instance " \
                    "WHERE %s_rating_type_id=? " \
                    "AND version=?" % tuple(3*[pdm_or_sp]) 
        db.execute(query, (type_id, version))
        row = db.cursor.fetchone()

        # Add rating instance if necessary
        if row is not None:
            id = row[0]
        else:
            query = "INSERT INTO %s_rating_instance " \
                        "(%s_rating_type_id, version, description) " \
                        "VALUES (?, ?, ?)" % tuple(2*[pdm_or_sp]) 
            db.execute(query, (type_id, version, description))
            query = "SELECT %s_rating_instance_id " \
                        "FROM %s_rating_instance " \
                        "WHERE %s_rating_type_id=? " \
                        "AND version=?" % tuple(3*[pdm_or_sp]) 
            db.execute(query,(type_id, version))
            id = db.cursor.fetchone()[0]
        db.close()

        # Return newly added rating instance id
        return id


def get_scaled_profile(profile, varprof):
    scaled = profile.copy()
    scaled /= np.sqrt(varprof)
    scaled -= scaled.mean()
    return scaled


def vonmises_coefficient(k,m):
    return scipy.special.ive(m,k)/scipy.special.ive(0,k)


def vonmises_values(k,mu,xs):
    distribution = scipy.stats.vonmises(k,scale=1./(2*np.pi))
    return distribution.pdf((xs-mu)%1)


def vonmises_histogram(k,mu,n,factor=2):
    if n % 2:
        raise ValueError("n (%d) must be even" % n)
    m = ((n*factor)//2+1)
    coeffs = vonmises_coefficient(k,np.arange(m)) * \
                np.exp(-2.0j*np.pi*mu*np.arange(m))*n*factor
    longhist = 1.0 + np.fft.irfft(coeffs*(np.exp(2.0j * \
                                   np.pi*1.0/n*np.arange(m))-1) / \
                   (np.maximum(np.arange(m), 1)*2.0j*np.pi))*n*factor
    return np.mean(np.reshape(longhist, (n,factor)), axis=-1)/factor

def multigaussfit_from_paramlist(params):
    comps = []
    for ii in range(1, len(params), 3):
        amp = params[ii]
        std = np.abs(params[ii+1])
        phs = params[ii+2]
        comps.append(dataproducts.MultiGaussComponent(amp, std, phs))
    fit = dataproducts.MultiGaussFit(offset=params[0], components=comps)
    return fit


def print_raters_list(verbosity=0):
    """Print the list of imported raters to stdout.
        
        Input:
            verbosity: If True, print description of raters.
                (Default: Don't be verbose.)

        Outputs:
            None
    """
    import textwrap
    import raters
    print "Number of raters registered: %d" % len(raters.registered_raters)
    for rater_name in raters.registered_raters:
        rater_module = getattr(raters, rater_name)
        rater = rater_module.Rater()
        print "'%s': %s (v%d)" % (rater_name, rater.long_name, rater.version)
        if verbosity:
            print ""
            for line in rater.description.split('\n'):
                print textwrap.fill(line, width=70)
            print "-"*25

def print_sp_raters_list(verbosity=0):
    """Print the list of imported single pulse raters to stdout.
        
        Input:
            verbosity: If True, print description of raters.
                (Default: Don't be verbose.)

        Outputs:
            None
    """
    import textwrap
    import sp_raters
    print "Number of single pulse raters registered: %d" % len(sp_raters.registered_raters)
    for rater_name in sp_raters.registered_raters:
        rater_module = getattr(sp_raters, rater_name)
        rater = rater_module.Rater()
        print "'%s': %s (v%d)" % (rater_name, rater.long_name, rater.version)
        if verbosity:
            print ""
            for line in rater.description.split('\n'):
                print textwrap.fill(line, width=70)
            print "-"*25
