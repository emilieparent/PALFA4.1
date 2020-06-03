import numpy as np
import psr_utils
from scipy.optimize import leastsq
from scipy import stats

import base
from rating_classes import pfd as pfdcls

class DMCurveRater(base.BaseRater):
    short_name = "dm_curve"
    long_name = "DM Curve"
    description = "Generate a simulated set of subbands using the 1D " \
                  "summed profile, dedisperse the real and simulated " \
                  "data at a set of trial DMs, then get a reduced " \
                  "chi-squared value for the difference between the two " \
                  "DM-vs-chisq curves. Pulsars with reasonable statistics " \
                  "(not crazy bright ones) should return a value near or " \
                  "below 1."
    version = 1

    rat_cls = pfdcls.PfdRatingClass()

    def _compute_rating(self, cand):
        """Return a rating for the candidate.
            Generate a simulated set of subbands using the 1D
            summed profile, dedisperse the real and simulated
            data at a set of trial DMs, then get a reduced
            chi-squared value for the difference between the two
            DM-vs-chisq curves. Pulsars with reasonable statistics
            (not crazy bright ones) should return a value near or
            below 1.
           

            Input:
                cand: A Candidate object to rate.

            Output:
                value: The rating value.
        """
        self.pfd = cand.get_from_cache('pfd')

        spec_index = self.find_spec_index()

        dms, dmcurve, fdmcurve = self.dm_curve_check(spec_index)

        return np.sum(self.dm_curve_diff(spec_index)**2) / (len(dms)-1)

    def dm_curve_check(self, spec_index=0.):
        # Sum the profiles in time
        profs = self.pfd.profs.sum(0)

        ### Generate simulated profiles ###

        # prof_avg: median profile value per subint per subband
        #  Sum over subint axis to get median per subband
        prof_avg = self.pfd.stats[:,:,4].sum(0)
        # prof_var: variance of profile per subint per subband
        #  Sum over subint axis to get variance per subband
        prof_var = self.pfd.stats[:,:,5].sum(0)
        # The standard deviation in each subband is proportional to the median
        # value of that subband.  Here we scale all subbands to equal levels.
        #scaled_vars = prof_var / prof_avg**2
        scaled_vars = np.array(np.abs(prof_var),dtype=np.float64)/prof_avg**2
	scaled_vars[scaled_vars<=0]=np.random.rand(1) ## This is new, to avoir errors when noise calculates <=0 inside the sqrt
        scaled_profs = (profs.T / prof_avg).T - 1.
        # The mean profile (after scaling) will be our "clean" profile--hardly
        # true in most cases, but we pretend it's noiseless for what follows
        scaled_mean_prof = scaled_profs.mean(0)
        # Extend this "clean" profile across the number of subbands
        sim_profs_clean = np.tile(scaled_mean_prof,\
            scaled_profs.shape[0]).reshape(scaled_profs.shape)
        # Scale these subbands according to the input spectral index
        spec_mult = (self.pfd.subfreqs/self.pfd.subfreqs[0])**spec_index
        spec_mult /= spec_mult.mean()
        sim_profs_spec = (sim_profs_clean.T*spec_mult).T
        # For consistency, set a seed for generating noise
        np.random.seed(1967)
        # Add white noise that matches the variance of the real subbands
        # on a subband by subband basis
        noise = np.random.normal(scale=np.sqrt(scaled_vars),size=scaled_profs.T.shape).T
        # sim_profs_noisy is the simulated equivalent of scaled_profs
        sim_profs_noisy = sim_profs_spec + noise
        # sim_profs_final is the simulated equivalent of profs
        sim_profs = ((sim_profs_noisy + 1.).T * prof_avg).T
        
        # The rest of this is essentially code from the prepfold.pfd class
        # in which we loop over DM values and see how strong a signal we
        # get by dedispersing at each of these values

        # Go to higher DMs than the original curve to try to better exclude noise
        DMs = np.linspace(self.pfd.dms[0], self.pfd.dms[-1]+4.*\
            (self.pfd.dms[-1]-self.pfd.dms[0]), len(self.pfd.dms))
        chis = np.zeros_like(DMs)
        sim_chis = np.zeros_like(DMs)
        subdelays_bins = self.pfd.subdelays_bins.copy()
        for ii, DM in enumerate(DMs):
            subdelays = psr_utils.delay_from_DM(DM, self.pfd.barysubfreqs)
            hifreqdelay = subdelays[-1]
            subdelays = subdelays - hifreqdelay
            delaybins = subdelays*self.pfd.binspersec - subdelays_bins
            new_subdelays_bins = np.floor(delaybins+0.5)
            for jj in range(self.pfd.nsub):
                profs[jj] = psr_utils.rotate(profs[jj], int(new_subdelays_bins[jj]))
                sim_profs[jj] = psr_utils.rotate(sim_profs[jj],\
                    int(new_subdelays_bins[jj]))
            subdelays_bins += new_subdelays_bins
            # The set of reduced chi2s like those in the prepfold plot
            # (should be the same if the same DMs are used)
            chis[ii] = self.pfd.calc_redchi2(prof=profs.sum(0), avg=self.pfd.avgprof)
            # The same thing but for our "simulated" data
            sim_chis[ii] = self.pfd.calc_redchi2(prof=sim_profs.sum(0), avg=self.pfd.avgprof)
        return DMs, chis, sim_chis

    def dm_curve_diff(self, spec_index):
        dms, dmcurve, fdmcurve = self.dm_curve_check(spec_index)
        # I'm not quite sure whether to use the variance of both curves or just
        # one of them...

        # This variance comes from basic error propagation through the
        # pfd.calc_redchi2 function
        var1 = 4./self.pfd.DOFcor * dmcurve

        #var2 = 4./self.pfd.DOFcor * fdmcurve 

        return (dmcurve - fdmcurve) / np.sqrt(var1)

    def find_spec_index(self):
        """
        Run a least squares fit to find the "spectral index" that gives the
        best match between real and simulated DM curves
        """
        fit = leastsq(self.dm_curve_diff, x0=0, maxfev=20)
        return fit[0][0]

Rater = DMCurveRater
