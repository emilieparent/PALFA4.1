"""
pipeline_utils.py

Defines utilities that will be re-used by multiple modules/scripts
in the pipeline package.

Emilie Parent, July 2019 (adapted from P. Lazarus) 
"""
import os
import os.path
import sys
import subprocess
import types
import traceback
import optparse
import time
import datetime
import string
import shutil
import glob

import debug
import jobtracker2
import config.basic 
import config.processing
import config.searching
import shutil
import tarfile
import astro_utils

import getpass
from pwd import getpwuid
from grp import getgrgid


def get_file_ownership(filename):
	return getpwuid(os.stat(filename).st_uid).pw_name

class PipelineError(Exception):
    """A generic exception to be thrown by the pipeline.
    """
    def __init__(self, *args, **kwargs):
        super(PipelineError, self).__init__(*args, **kwargs)
        exctype, excval, exctb = sys.exc_info()
        if (exctype is not None) and (excval is not None) and \
                (exctb is not None):
            self.orig_exc_info = exctype, excval, exctb

    def __str__(self):
        msg = super(PipelineError, self).__str__()
        if 'orig_exc_info' in self.__dict__.keys():
            msg += "\n\n========== Original Traceback ==========\n"
            msg += "".join(traceback.format_exception(*self.orig_exc_info))
            msg += "\n(See PipelineError traceback above)\n"
        if msg.count("\n") > 100:
            msg = string.join(msg.split("\n")[:50],"\n")
        return msg


def get_fns_for_jobid(jobid):
    """Given a job ID number, return a list of that job's data files.

        Input:
            jobid: The ID number from the job-tracker DB to get files for.
        
        Output:
            fns: A list of data files associated with the job ID.
    """

    query = "SELECT filename " \
            "FROM files, job_files " \
            "WHERE job_files.file_id=files.id " \
                "AND job_files.job_id=%d" % jobid
    rows = jobtracker2.query(query)
    fns = [str(row['filename']) for row in rows]
    return fns


def clean_up(jobid):
    """Deletes raw files for a given job ID.

        Input:
            jobid: The ID corresponding to a row from the job_submits table.
                The files associated to this job will be removed.
        Outputs:
            None
    """
    fns = get_fns_for_jobid(jobid)
    for fn in fns:
        remove_file(fn)


def archive_logs():
    """
	New - Emilie Parent, July 2019    
	    Removes the logs files for a jobid that was sucessfully 
            uploaded  to the database at Cornell. 
            This is to avoid accumulating files, since the number of 
            files limit is low on Beluga.
	Input:
            jobid: The ID corresponding to a row from the job_submits table.
                   The files associated to this job will be removed. (job_id,
	           unlike id, uniquely defines a beam.
	Output:
            logs_to_del: List of logs.ER and logs.OU files that can be deleted.
 
    """
    curdir = os.getcwd()
    os.chdir(config.basic.qsublog_dir)

    tar = tarfile.open('/project/ctb-vkaspi/PALFA4/archived_logs2.tar','a')

    query = "SELECT * FROM job_submits WHERE status IN "\
	    "('uploaded','processed') AND updated_at>'2019-08-01 00:00:00'"
    rows = jobtracker2.query(query)
    queue_ids = [str(rows[i]['queue_id']+'.*') for i in range(len(rows))]

    for q in queue_ids:
	f = glob.glob(q)
	if len(f)==2:
		tar.add(f[0])
		tar.add(f[1])
		print "Added a job's logs to the archived logs"
		os.remove(f[0])
		os.remove(f[1])

    tar.close()
    os.chdir(curdir)


def remove_logs(date='2019-07-01'):
    """
	New - Emilie Parent, July 2019    
	    Removes the logs files for a jobid that was sucessfully 
            uploaded  to the database at Cornell. 
            This is to avoid accumulating files, since the number of 
            files limit is low on Beluga.
	Input:
            date: look in the database for logs of jobs updated after that date (default = July 2019, i.e. on Beluga) 
	Output:
	    -
 
    """
    curdir = os.getcwd()
    os.chdir(config.basic.qsublog_dir)

    query1 = "SELECT id FROM jobs WHERE status='processed' AND details in ('Ready for upload','Processed without errors') and updated_at>'%s 00:00:00' ORDER BY id DESC"%date
    query2 = "SELECT id FROM jobs WHERE status='uploaded' and updated_at>'%s 00:00:00' ORDER BY id DESC"%date
    rows = jobtracker2.query(query1)
    rows +=jobtracker2.query(query2)
    job_ids = [str(r[0]) for r in rows]
    k=0
    for j in job_ids:
	jsub = jobtracker2.query("SELECT * FROM job_submits where job_id=%s ORDER BY id DESC"%j)
	queue_ids = [str(s['queue_id'])+'.*' for s in jsub]
	for q in queue_ids:
		f = glob.glob(q)
		if len(f)==2:
			os.remove(f[0])
			os.remove(f[1])
			k+=1
		elif len(f)==1:
			os.remove(f[0])
			k+=1
    print "Removed logs of %d submitted jobs"%k
    os.chdir(curdir)

def change_job_status_jtdb(ROW, new_status, details='Errors with results files'):
    """
	New - Emilie Parent, July 2019
	Called by move_results when the following error is encountered:
		a job has the status "processed" in the jobtracking db, 
		but it's result directory in scratch is empty. Weird. 
	Input:
		ROW:        a row from job_submits table
		new_status: 'failed' should be provided as new status when used by move_results()
		id:         id of the job (not job_id or queue_id)

    """
    
    queries = []
    
    if new_status=='processing_failed' or new_status=='submission_failed' or new_status=='precheck_failed' or new_status=='upload_failed':
	jobs_status='failed'
    elif new_status=='processed' or 'uploaded':
	jobs_status=new_status
    elif new_status=='running':
	jobs_status='submitted'

    queries.append("UPDATE jobs " \
             "SET status='%s', updated_at='%s', details='%s' " \
             "WHERE id=%d" % (jobs_status, jobtracker2.nowstr(),details, ROW['job_id']))

    queries.append("UPDATE job_submits " \
             "SET status='%s', updated_at='%s', details='%s' " \
             "WHERE id=%d" % (new_status, jobtracker2.nowstr(), details, ROW['id']))
    jobtracker2.query(queries)


def go_copy_results(row, end_name, new_dir):
    """
	Step 2: copy the results from scratch to projects, and then remove initial copy on scratch
    	This should solve the problem of having the username as the name of the group on projects 
    """
    allowed_usr = config.processing.users
    current_usr = getpass.getuser()

    try:
	cp_cmd = "cp -r "+str(row['output_dir'])+'/* '+new_dir+'/'
	subprocess.call(cp_cmd,shell=True)
	rm_cmd = "rm -r "+str(row['output_dir'])
	subprocess.call(rm_cmd,shell=True)

	jobtracker2.query("update job_submits set output_dir='%s' where job_id='%s'"%(new_dir,str(row['job_id'])))
	change_job_status_jtdb(ROW=row, new_status = 'processed', details='Ready for upload')
	print "Moved results to project area, and updated jtdb for: ", row['output_dir'].split('/results/')[-1]

	path1 = os.path.split(new_dir)[0]
	path2 = os.path.split(path1)[0]
	path3 = os.path.split(path2)[0]

	owner_newdir = get_file_ownership(new_dir)
	owner1 = get_file_ownership(path1)
	owner2 = get_file_ownership(path2)
	owner3 = get_file_ownership(path3)

	paths  = [new_dir,     path1,  path2,  path3]
	owners = [owner_newdir,owner1, owner2, owner3]

	for i in range(len(paths)):
		if owners[i]==current_usr:
			if i==0:
				for user in allowed_usr:
					cmd = "setfacl -R -m u:%s:rwx %s"%(user,paths[i])
					subprocess.call(cmd,shell=True)
			else:
				for user in allowed_usr:
					cmd = "setfacl -m u:%s:rwx %s"%(user,paths[i])
					subprocess.call(cmd,shell=True)	

	return '1'

    except: 
	print "Error while moving results from scratch to project ",\
		"for beam: ", end_name, "queue_id=",row['queue_id']
	path, dirs, files = next(os.walk(row['output_dir']))
	if len(files)==0:
		print "\t","Job's status is 'processed', scratch result directory exists but is empty.."
		print "\t",'Changing job status in jobtracking database for "failed". To be retried.'
		change_job_status_jtdb(ROW=row,new_status='failed')
		return '0'
				
	else:
		print "\t","Unknown error.. check it out manually."
		print "Exiting "
		sys.exit()




def move_results():
    """ 
	New: Moves all result directories for which the job was successfully processed
	from config.processing.base_results_directory to config.processing.base_final_results_dir 

	 by Emilie Parent, July 2019   
    """

    query = "SELECT * FROM job_submits WHERE status='processed' AND output_dir like '/scratch%'"
    rows = jobtracker2.query(query)
    
    current_usr = getpass.getuser()
    j = 0
    skipped = 0
    N = len(rows)
    if rows:
	print "Will attempt to move %d job results"%len(rows)
    	for i,row in enumerate(rows):
		owner_results = get_file_ownership(row['output_dir'])
		if owner_results != current_usr:
			print "	Wrong user: Owner is %s, current user is %s. Skipping."%(owner_results,current_usr)
			skipped+=1
			continue
		end_name = str(row['output_dir']).split('results/')[1]
		while end_name.endswith('/'):
			end_name = end_name[:-1]
 
		new_dir = str(config.processing.base_final_results_dir)+end_name
		if i%10==0:
			print i,' / ',N
		# The results have already been moved to projects area
		if os.path.exists(row['output_dir']) and ('projects' in row['output_dir']):
			pass

		# there is an existing result directory in the scratch 
		# and the job is processed : to be moved to project area
		elif os.path.exists(row['output_dir']) and ('scratch' in row['output_dir']):
			# Step 1: try to create the new result directory in projects area
			try:
				tar = glob.glob(row['output_dir']+'/p2030*00.tgz')
				if len(tar)==0 and not os.path.exists(new_dir):
					#Problem: did not produce all files needed
					print "Incomplete processing for %s, " \
					"changing status to 'failed'"%end_name 
					change_job_status_jtdb(row, new_status='processing_failed', details='Errors with results files')
					pass
					
				os.makedirs(new_dir)
				add = go_copy_results(row, end_name, new_dir)
				j+=int(add)

			except:
				print "Beam ",end_name,"  already has a result ",\
					"directory in project area? Checking .. "
			
				if not os.path.exists(new_dir):
					print "\t","No, it does not exist.",\
						"Check manually for possible problems,",\
					"queue_id = ",row['queue_id']
					# that's a weird case: something is wrong
					tar1 = glob.glob(row['output_dir']+'/p2030*00.tgz')
					tar2 = glob.glob(row['output_dir']+'/p2030*.fits')
					if (len(tar1)==0 or len(tar2)==0) and not os.path.exists(new_dir):
						print "Incomplete processing for %s, " \
						"changing status to 'failed'"%end_name
						change_job_status_jtdb(row, new_status='processing_failed', details='Errors with results files')
						
					

				else:
					path, dirs, files = next(os.walk(new_dir))
					if len(files)==0:
						print "\t","Directory already exists in projects ",\
						"but is empty: copying results"
						add = go_copy_results(row, end_name, new_dir)
						j+=int(add)
					else:
						# should check why there is stuff in that directory, 
						# and yet the original on scratch has not been deleted.. 
						print "\t","Directory already exists AND not empty: ",\
							"not copying over beam ",end_name

		# there is no result directory in either projects or scratch
		# already uploaded, or some problem happended.
		else:
			print "Beam ",end_name,": no existing result directory in /scratch"
			change_job_status_jtdb(row, new_status='processing_failed', details='Errors with results files')
			pass
			
    if j>0 or skipped>0:
    	print '\n',"Moved %s processed beams to project area (skipped %s). "%(str(j),str(skipped)),"\t -> Ready for upload.",'\n'



def remove_file(fn):
    """Delete a file (if it exists) and mark it as deleted in the 
        job-tracker DB.

        Input:
            fn: The name of the file to remove.

        Outputs:
            None
    """
    if os.path.exists(fn):
        os.remove(fn)
        print "Deleted: %s" % fn
    jobtracker2.query("UPDATE files " \
                     "SET status='deleted', " \
                         "updated_at='%s', " \
                         "details='File was deleted' " \
                     "WHERE filename='%s'" % \
                     (jobtracker2.nowstr(), fn))


def can_add_file(fn, verbose=False):
    """Checks a file to see if it should be added to the 'files'
        table in the jobtracker DB.

        Input:
            fn: The file to check.
            verbose: Print messages to stdout. (Default: be silent).

        Outputs:
            can_add: Boolean value. True if the file should be added. 
                    False otherwise.
    """
    import datafile
    try:
        datafile_type = datafile.get_datafile_type([fn])
    except datafile.DataFileError, e:
        if verbose:
            print "Unrecognized data file type: %s" % fn
        return False
    parsedfn = datafile_type.fnmatch(fn)
    if parsedfn.groupdict().setdefault('beam', '-1') == '7':
        if verbose:
            print "Ignoring beam 7 data: %s" % fn
        return False
    # Check if file is already in the job-tracker DB
    files = jobtracker2.query("SELECT * FROM files " \
                             "WHERE filename LIKE '%%%s'" % os.path.split(fn)[-1])
    if len(files):
        if verbose:
            print "File is already being tracked: %s" % fn
        return False

    # Check if file has a corresponding custom zaplist
    if not config.processing.use_default_zaplists \
       and not find_zaplist_in_tarball(fn,verbose=verbose):
        return False

    return True


def execute(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, dir=None): 
    """Execute the command 'cmd' after logging the command
        to STDOUT. Execute the command in the directory 'dir',
        which defaults to the current directory is not provided.

        Output standard output to 'stdout' and standard
        error to 'stderr'. Both are strings containing filenames.
        If values are None, the out/err streams are not recorded.
        By default stdout is subprocess.PIPE and stderr is sent 
        to sys.stderr.

        Returns (stdoutdata, stderrdata). These will both be None, 
        unless subprocess.PIPE is provided.
    """
    # Log command to stdout
    if debug.SYSCALLS:
        sys.stdout.write("\n'"+cmd+"'\n")
        sys.stdout.flush()

    stdoutfile = False
    stderrfile = False
    if type(stdout) == types.StringType:
        stdout = open(stdout, 'w')
        stdoutfile = True
    if type(stderr) == types.StringType:
        stderr = open(stderr, 'w')
        stderrfile = True
    
    # Run (and time) the command. Check for errors.
    pipe = subprocess.Popen(cmd, shell=True, cwd=dir, \
                            stdout=stdout, stderr=stderr)
    (stdoutdata, stderrdata) = pipe.communicate()
    retcode = pipe.returncode 
    if retcode < 0:
        raise PipelineError("Execution of command (%s) terminated by signal (%s)!" % \
                                (cmd, -retcode))
    elif retcode > 0:
        raise PipelineError("Execution of command (%s) failed with status (%s)!" % \
                                (cmd, retcode))
    else:
        # Exit code is 0, which is "Success". Do nothing.
        pass
    
    # Close file objects, if any
    if stdoutfile:
        stdout.close()
    if stderrfile:
        stderr.close()

    return (stdoutdata, stderrdata)

def get_modtime(file, local=False):
    """Get modification time of a file.

        Inputs:
            file: The file to get modification time for.
            local: Boolean value. If true return modtime with respect 
                to local timezone. Otherwise return modtime with respect
                to GMT. (Default: GMT).

        Outputs:
            modtime: A datetime.datetime object that encodes the modification
                time of 'file'.
    """
    if local:
        modtime = datetime.datetime(*time.localtime(os.path.getmtime(file))[:6])
    else:
        modtime = datetime.datetime(*time.gmtime(os.path.getmtime(file))[:6])
    return modtime


def get_zaplist_tarball(force_download=False, no_check=False, verbose=False):
    """Download zaplist tarball. If the local version has a
        modification time equal to, or later than, the version
        on the FTP server don't download unless 'force_download'
        is True.

        Input:
            force_download: Download zaplist tarball regardless
                of modification times.
            no_check: Do not check for updated tarball. Just
                return name.
            verbose: If True, print messages to stdout.
                (Default: Be silent)

        Outputs:
            Name with full path of zaplist tarball.
    """
    import CornellFTP


    if config.searching.use_radar_clipping:
        zaptar_basename = "PALFA4_zaplists_noradar.tar.gz"
    else:
        zaptar_basename = "zaplists.tar.gz"
    
    zaptarfile = os.path.join(config.processing.zaplistdir, zaptar_basename)

    if no_check:
        return zaptarfile

    cftp = CornellFTP.CornellFTP()

    ftpzappath = "/zaplists/" + zaptar_basename
    getzap = False
    if force_download:
        if verbose:
            print "Forcing download of zaplist tarball"
        getzap = True
    elif not os.path.exists(zaptarfile):
        if verbose:
            print "Zaplist tarball doesn't exist, will download"
        getzap = True
    elif cftp.get_modtime(ftpzappath) > get_modtime(zaptarfile):
        if verbose:
            print "Zaplist on FTP server is newer than local copy, will download"
        getzap = True

    cftp.close()

    if getzap:
        list_fn = zaptar_basename.replace('.tar.gz','.list')
        temp_list_fn = zaptar_basename.replace('.tar.gz','_dl.list')
        zaplistdir = config.processing.zaplistdir

        temp_zaplistfn = os.path.join(zaplistdir,'zaplists_dl.tar.gz')

        # Download the file from the FTP
        CornellFTP.pget(ftpzappath, temp_zaplistfn)

        # Make text list of zaplist tarball contents to speed up
        # finding of zaplists in tarball
        zaptar = tarfile.open(temp_zaplistfn, mode='r')
        names = zaptar.getnames()

        zaplistf = open(os.path.join(zaplistdir, \
                        temp_list_fn),'w')
        for name in names:
            zaplistf.write(name+'\n')
        
        zaplistf.close()
        zaptar.close()

        os.rename(temp_zaplistfn, zaptarfile)
        os.rename(os.path.join(zaplistdir,temp_list_fn), \
                  os.path.join(zaplistdir,list_fn)) 
        
    else:
        # Do nothing
        pass

    return zaptarfile

def find_zaplist_in_tarball(filename, verbose=False):
    """Find the name of the zaplist for a given raw data filename.
        Searches the 'zaplists_tarball.list' textfile for the name
        of the zaplist corresponding to the raw data file.        

        Input: filename - name of the raw data file.
 
        Output: zaplist - name of the zaplist in the tarball.
    """
    import datafile
    fns = [ filename ]
    filetype = datafile.get_datafile_type(fns)
    parsed = filetype.fnmatch(fns[0]).groupdict()
    if 'date' not in parsed.keys():
        parsed['date'] = "%04d%02d%02d" % \
                            astro_utils.calendar.MJD_to_date(int(parsed['mjd']))

    zaplist_tarball_fn = get_zaplist_tarball(no_check=True)

    if verbose:
        print "Looking for zaplist for %s in %s..." % \
               ( filename, os.path.basename(zaplist_tarball_fn) )

    customzapfns = []
    # First, try to find a custom zaplist for this specific data file
    customzapfns.append("%s.%s.%s.b%s.%s.zaplist" % \
                        (parsed['projid'], parsed['date'], parsed['source'], \
                         parsed['beam'], parsed['scan']))
    # Next, try to find custom zaplist for this beam
    customzapfns.append("%s.%s.b%s.zaplist" % \
                        (parsed['projid'], parsed['date'], parsed['beam']))
    # Try to find custom zaplist for this MJD
    customzapfns.append("%s.%s.all.zaplist" % (parsed['projid'], parsed['date']))

    zaplistf = open(zaplist_tarball_fn.replace('.tar.gz','.list'),'r')
    names = zaplistf.readlines()
    zaplistf.close()

    for customzapfn in customzapfns:
        matches = [name for name in names \
                    if name.endswith(customzapfn+'\n')]
        if matches:
            zaplist = matches[0].rstrip('\n')
            if verbose:
                print "Found zaplist",zaplist
            return zaplist
        else:
            # The member we searched for doesn't exist, try next one
            pass
    else:
        # No custom zaplist found.
        if verbose:
            print "No zaplist found."
        return None

class PipelineOptions(optparse.OptionParser):
    def __init__(self, *args, **kwargs):
        optparse.OptionParser.__init__(self, *args, **kwargs)
       
    def parse_args(self, *args, **kwargs):
        # Add debug group just before parsing so it is the last set of
        # options displayed in help text
        self.add_debug_group()
        return optparse.OptionParser.parse_args(self, *args, **kwargs)

    def add_debug_group(self):
        group = optparse.OptionGroup(self, "Debug Options", \
                    "The following options turn on various debugging " \
                    "features in the pipeline. Multiple debugging " \
                    "options can be provided.")
        group.add_option('-d', '--debug', action='callback', \
                          callback=self.debugall_callback, \
                          help="Turn on all debugging modes. (Same as --debug-all).")
        group.add_option('--debug-all', action='callback', \
                          callback=self.debugall_callback, \
                          help="Turn on all debugging modes. (Same as -d/--debug).")
        for m, desc in debug.modes:
            group.add_option('--debug-%s' % m.lower(), action='callback', \
                              callback=self.debug_callback, \
                              callback_args=(m,), \
                              help=desc)
        self.add_option_group(group)

    def debug_callback(self, option, opt_str, value, parser, mode):
        debug.set_mode_on(mode)

    def debugall_callback(self, option, opt_str, value, parser):
        debug.set_allmodes_on()
