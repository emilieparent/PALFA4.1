#!/usr/bin/env python
import sys
import re
import os
import copy
import numpy as Num
import subprocess 

from presto import candidate_sigma
import sifting

"""
Slightly modified version of sifting.py (PRESTO). 

"""
# Longest period candidates to consider (s)
long_period = 30.0
# Shortest period candidates to consider (s)
short_period = 0.05



fund_re = re.compile("^\d")
harms_re = re.compile("^[ ]\d")
DM_re = re.compile("DM(\d+\.\d{2})")
prelim_reject = False


def cmp_snr(self, other):
    retval = -cmp(self.snr, other.snr)
    return retval
#==========================================

class FFACandidate(sifting.Candidate):
    def __init__(self,candnum, p, snr, dt ,binn, dm, DMstr, filename, T, final_cands=False):
        self.path, self.filename = os.path.split(filename)
	# to identify the filename, replace "_precands.ffa" by ".dat"
	self.filename= self.filename.split('_precands.ffa',1)[0]
	self.filename = self.filename + '.dat'
	if final_cands :
		#if final_cands, the candnum is replaced by the filename
		self.candnum = (candnum)
	else:	
		self.candnum = int(candnum)
	self.p = p
	#for the ffa, snr and sigma are the same thing.
        self.snr = snr
	self.sigma = snr
        self.numharm = 1  
	self.dt = dt
        self.f = 1.0/p
        self.T = T
	self.r = binn
        self.DMstr = DMstr
        self.DM = float(dm)
	self.hits = []
        self.note = ""
    def add_as_hit(self, other):
        self.hits.extend([(other.DM,other.snr,other.sigma)])

    def __str__(self):
        cand = self.filename + '    ' + `self.candnum`
        return "%-65s   %3.1f  %5.2f   %5.1f   %5.1f"% (cand, self.p*1000, self.snr, self.DM, self.dt*1000)

#==========================================

class FFACandlist(sifting.Candlist):
    def __init__(self, cands=None, trackbad=False, trackdupes=False):
        if cands is None:
            self.cands = []
        else:
            self.cands = cands
        self.trackbad = trackbad # Should we keep track of bad candidates
        self.trackdupes = trackdupes # Should we keep track of duplicates
        # Set default badlists
# REMOVED HARMONIC RELATED STUFF
        self.badlists = {'knownbirds': [], \
                         'longperiod': [], \
                         'shortperiod': [], \
                         'threshold': [], \
                         'dmproblem': []}
        self.duplicates = []

    def plot_ffa_summary(self, usefreqs=True):
        """Produce a plot summarizing the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Get all candidates and sort by sigma
        allcands = self.get_all_cands()
        sigmas = Num.array([c.sigma for c in allcands])
        isort = sigmas.argsort()
        sigmas = Num.array(sigmas[isort])

        if usefreqs:
            xdata = Num.array([c.f for c in allcands])[isort]
            xlabel = "Freq (Hz)"
            xscale = "log"
        else:
            xdata = Num.array([c.p for c in allcands])[isort]
            xlabel = "Period (s)"
            xscale = "loglin"

        dms = Num.array([c.DM for c in allcands])[isort]
        numharms = Num.array([c.numharm for c in allcands])[isort]

        # Plot the all candidates 
        scatt = plt.scatter(xdata, dms, s=sifting.sigma_to_size(sigmas), \
                                c=Num.log2(numharms), \
                                marker='o', alpha=0.7, zorder=-1) 
        plt.set_cmap("Spectral") 
  
        # Add colorbar 
        #fmtr = matplotlib.ticker.FuncFormatter(lambda x, pos: "%d" % 2**x)
        #cax = plt.axes((0.18, 0.06, 0.67, 0.035))
        #cb = plt.colorbar(scatt, cax=cax, ticks=(0,1,2,3,4), format=fmtr, \
        #                    orientation="horizontal")
        #cb.set_label("Num harmonics summed") 
        
        plt.axes(ax) # Set scatter plot's axes as current
        plt.xscale(xscale)
        plt.xlabel(xlabel)
        mindm = Num.min(dms)
        maxdm = Num.max(dms)
        dmrange = Num.ptp(dms)

        # Use log-scale y-axis if max DM > 2000
        yscale = "log" if maxdm > 2000.0 else "linear"
        plt.yscale(yscale)

        if yscale is "log":
            plt.ylim(1.0, maxdm+0.1*dmrange)
        else:
            plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)

        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        if not usefreqs:
            plt.gca().xaxis.set_ticks(Num.concatenate((\
                                        Num.logspace(-4,0,4, endpoint=False), \
                                        Num.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(Num.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(xdata)/5.0), \
                        min(long_period+0.5, max(xdata)+0.5))
        else:
            plt.xlim(min(xdata)/5.0, max(xdata)*2.0)
            
        ax.format_coord = lambda x,y: "x=%g, y=%g" % (x,y)
        return fig
    
    def plot_ffa_rejects(self, usefreqs=True):
        """Produce a plot showing why candidates were rejected by
            the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Plot bad candidates
        candlists = [self.badlists['knownbirds'], self.badlists['longperiod'], \
                     self.badlists['shortperiod'], self.badlists['threshold'], \
                     self.badlists['dmproblem'], \
                     self.cands, self.duplicates]
        labels = ['Known birdires', 'Long period', 'Short period', \
                    'Threshold', 'DM problem', 'Good cands', 'Hits']
        colours = ['#FF0000', '#800000', '#008000', '#00FF00', \
                    '#800080', 'r', 'k']
        markers = ['o', 'o', 'o', 'o', 'o', 'x', 's']
        zorders = [-2, -2, -2, -2, -2, 0, 0]
        sizes = [50, 50, 50, 50, 50, 100, 10]
        fixedsizes = [0, 0, 0, 0, 0, 1, 1]
        lws = [1,1,1,1,1,2,1]
        handles = []
        for cands, colour, marker, zorder, size, fixedsize, lw in \
                zip(candlists, colours, markers, zorders, sizes, fixedsizes, lws):
            if len(cands):
                sigmas = Num.array([c.sigma for c in cands])
                isort = sigmas.argsort()
                sigmas = sigmas[isort]
                if usefreqs:
                    xdata = Num.array([c.f for c in cands])[isort]
                    xlabel = "Freq (Hz)"
                    xscale = "log"
                else:
                    xdata = Num.array([c.p for c in cands])[isort]
                    xlabel = "Period (s)"
                    xscale = "loglin"
                dms = Num.array([c.DM for c in cands])[isort]
                
                # Plot the candidates
                if fixedsize:
                    plt.scatter(xdata, dms, s=size, lw=lw, \
                                c=colour, marker=marker, alpha=0.7, zorder=zorder)
                else:
                    plt.scatter(xdata, dms, s=sigma_to_size(sigmas), lw=lw, \
                                c=colour, marker=marker, alpha=0.7, zorder=zorder)
            handles.append(plt.scatter([0], [0], s=size, c=colour, \
                                    marker=marker, alpha=0.7))

        fig.legend(handles, labels, 'lower center', \
                        prop={'size':'x-small'}, ncol=4)

        plt.xscale(xscale) 
        plt.xlabel(xlabel)

        alldms = Num.array([c.DM for c in self.get_all_cands()])
        mindm = Num.min(alldms)
        maxdm = Num.max(alldms)
        dmrange = Num.ptp(alldms)

        # Use log-scale y-axis if max DM > 2000
        yscale = "log" if maxdm > 2000.0 else "linear"
        plt.yscale(yscale)

        if yscale is "log":
            plt.ylim(1.0, maxdm+0.1*dmrange)
        else:
            plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)

        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        if not usefreqs:
            all_xdata = Num.array([c.p for c in self.get_all_cands()])
            plt.gca().xaxis.set_ticks(Num.concatenate((\
                                        Num.logspace(-4,0,4, endpoint=False), \
                                        Num.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(Num.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(all_xdata)/5.0), \
                        min(long_period+0.5, max(all_xdata)+0.5))
        else:
            all_xdata = Num.array([c.f for c in self.get_all_cands()])
            plt.xlim(min(all_xdata)/5.0, max(all_xdata)*2.0)

        return fig
    #def plot_goodcands(self, usefreqs=True):
	#print "plot_goodcands: Not available for FFA candidates"
    #def reject_threshold(self, sigma_threshold=None,c_pow_threshold=None):
	#print "reject_threshold: Not available for FFA candidates"
    #def reject_harmpowcutoff(self, harm_pow_cutoff=None):
	#print "reject_harmpowcutoff: Not available for FFA candidates"
    #def reject_rogueharmpow(self):
	#print "reject_rogueharmpow: Not available for FFA candidates"
    #def print_cand_summary(self, summaryfilenm=None):
	#print "print_cand_summary: Not available for FFA candidates, should be thought."
    #def write_cand_report(self, reportfilenm=None):
	#print "write_cand_report: Not available for FFA candidates, should be thought."


    def remove_duplicate_candidates(self, verbosity=1):
        """Remove lower-significance 'duplicate' (i.e. same period)
            candidates from a list of candidates.  For the highest
            significance candidate, include a list of the DMs (and SNRs)
            of all the other detections.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        if verbosity >= 1:
            print "  Sorting the %d candidates by frequency..." % \
                        self.get_numcands()
        self.cands.sort(sifting.cmp_freq)
        if verbosity >= 1:
            print "  Searching for dupes..."
        ii = 0
        # Find any match
        while ii < self.get_numcands():
            jj = ii + 1
            if jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < sifting.r_err:
                # Find others that match
                jj += 1
                while jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < sifting.r_err:
                    jj += 1
                matches = self.cands[ii:jj]
                matches.sort(cmp_snr)
                bestindex = self.cands.index(matches[0])
                # flag the duplicates
                bestcand = self.cands[bestindex]
                # Add other matching cands as hit of highest-sigma cand
                for matchind in reversed(range(ii, jj)):
                    if matchind == bestindex:
                        # The current candidate is the highest-sigma cand, keep it
                        continue
                    match = self.cands[matchind]
                    bestcand.add_as_hit(match)
                    match.note = "This candidate is a duplicate of %s:%d" % \
                                (bestcand.filename, bestcand.candnum)
                    self.mark_as_duplicate(matchind)
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (match.filename, match.candnum, matchind)
                        print "    %s" % match.note
                # If the best candidate isn't at the same freq
                # as ii, then it's possible even more hits should
                # be added. So we don't increment the index
                # (note that the best cand has moved into position ii).
            else:
                ii += 1 # No candidates to be added as hits, move on
        if verbosity >= 1:
            print "Found %d candidates." % self.get_numcands()
        self.cands.sort(cmp_snr)
     

    def remove_harmonics(self, verbosity=1):
 	"""Remove the candidates that are lower significance harmonics
            of other candidates from the candlist.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Note:  should probably put the harmonics into the fundamental as hits (use sets)
        numremoved = 0
        self.cands.sort(cmp_snr)
        f_err = sifting.r_err/self.cands[0].T
        if verbosity >= 1:
            print "\nSearching for duplicate harmonics..."
        ii = 0
        while 1:
            fundcand = self.cands[ii]
	    if len(self.cands) <= 2:
		print "Too few candidates to sift. Skipping" 
		break
            jj = len(self.cands) - 1
            zapj = 0
            while 1:
                harmcand = self.cands[jj]
                for factor in Num.arange(1.0, 17.0):
                    if Num.fabs(fundcand.f - harmcand.f*factor) < f_err*factor:
                        zapj = 1
                        harmstr = "1/%dth" % factor
                    elif Num.fabs(fundcand.f - harmcand.f/factor) < f_err/factor:
                        zapj = 1
                        if factor==2.0:
                            harmstr = "%dnd" % factor
                        else:
                            harmstr = "%dth" % factor
                    if zapj:
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%s) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        harmstr, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        break
                # Check a few other common ratios
                for numer,denom in zip([3.0, 5.0, 2.0, 4.0, 5.0, \
                                        3.0, 5.0, 2.0, 3.0, 4.0],
                                       [2.0, 2.0, 3.0, 3.0, 3.0, \
                                        4.0, 4.0, 5.0, 5.0, 5.0]):
                    factor = numer/denom
                    if Num.fabs(fundcand.f-harmcand.f*factor) < f_err*factor:
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%d/%dth) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        denom, \
                                        numer, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        harmstr = "%d/%dth" % (denom, numer)
                        zapj = 1
                        break
                if zapj:
                    harmcand.note = "This candidate (P=%.4f s, DM=%.2f) is " \
                                    "a harmonic (%s) of %s:%d " \
                                    "(P=%.4f s, DM=%.2f)." % \
                                (harmcand.p, harmcand.DM, harmstr, \
                                    fundcand.filename, fundcand.candnum, \
                                    fundcand.p, fundcand.DM)
                    numremoved += 1
                    self.mark_as_bad(jj, 'harmonic')
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (harmcand.filename, harmcand.candnum, jj)
                        print "    %s" % harmcand.note
                    zapj = 0
                jj -= 1
                if jj == ii:
                    break
            ii += 1
            if ii >= len(self.cands) - 1:
                break
        if verbosity >= 1:
            print "Removed a total of %d harmonics.\n" % numremoved


    def remove_DM_problems(self, numdms, dmlist, low_DM_cutoff, verbosity=1):
 	"""Remove the candidates where any of the following are true:
            1) The number of hits is < numdms
            2) The highest S/N candidate occurs below a DM of low_DM_cutoff
            3) The minimum difference in DM indices between the hits is > 1

            Inputs:
                numdms: The minimum number of hits for a good candidate.
                dmlist: List of DMs.
                low_DM_cutoff: The lowest DM possible for a good candidate.
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Create a dictionary where the key is the dmstr 
        # and the values are the index
        dmdict = {}
        dms = Num.unique([float(dm) for dm in dmlist])
        dmstrs = ['%.2f'%dm for dm in dms]
        dmdict = dict(zip(dmstrs, range(len(dms))))
        numremoved = 0
        num_toofew = 0
        num_toolow = 0
        num_gaps = 0
        self.cands.sort(cmp_snr)
        for ii in reversed(range(len(self.cands))):
            currcand = self.cands[ii]

            # Remove all the candidates with DM < low_DM_cutoff
            if float(currcand.DM) <= low_DM_cutoff:
                numremoved += 1
                num_toolow += 1
                #currcand.note = "Hit with max sigma (%g) has dm (%.2f) " \
                #                "<= low DM cutoff (%.2f) " % \
                #                    (hitsigma, hitdm, low_DM_cutoff)
                self.mark_as_bad(ii, 'low_dmproblem')
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            # Remove all the candidates without enough DM hits
            if len(currcand.hits) < numdms:
                numremoved += 1
                num_toofew += 1
                currcand.note = "Candidate has only %d DM hits. This is less " \
                                "than minimum for 'good' cands (%d hits)" % \
                                (len(currcand.hits), numdms)
                self.mark_as_bad(ii, 'numdms_dmproblem')
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            imax = Num.argmax(Num.array([hit[2] for hit in currcand.hits]))
            hitdm, hitsnr, hitsigma = currcand.hits[imax]
	    

            # Remove all the candidates where there are no hits at consecutive DMs
            if len(currcand.hits) > 1:
                currcand.hits.sort(sifting.cmp_dms)
                dm_indices = Num.asarray([dmdict["%.2f"%currcand.hits[jj][0]] for jj in range(len(currcand.hits))])
                min_dmind_diff = min(dm_indices[1:] - dm_indices[:-1])
                if min_dmind_diff > 1:
                    numremoved += 1
                    num_gaps += 1
                    currcand.note = "DM list of hits has gaps (i.e. " \
                                    "consecutive DMs don't have hits)."
                    self.mark_as_bad(ii, 'no_consecutive_dmproblem')
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (currcand.filename, currcand.candnum, ii)
                        print "    %s" % currcand.note
                    continue

        if verbosity >= 1:
            print "Removed %d candidates with DM problems.\n" % numremoved
        if verbosity >= 2:
            print "  # with too few hits:", num_toofew
            print "  # with peak SNR too low:", num_toolow
            print "  # with gaps in DM hits:", num_gaps
        
    def reject_knownbirds(self, known_birds_f):
        """Find and remove candidates conincident with known birds.

            Inputs:
                known_birds_f: A list of tuples containing bad frequencies
                    and widths. The tuples should contain
                        (<bad freq (Hz)>, <one-sided width (Hz)>)
                    (Default: Globally defined "known_birds_f")
                known_birds_p: A list of tuples containing bad peridocities
                    and widths. The tuples should contain
                        (<bad freq (ms)>, <one-sided width (ms)>)
                    (Default: Globally defined "known_birds_p")

            Outputs:
                None
        """
	known_birds_p = (1./known_birds_f[0],1./known_birds_f[1])
	bird_removed = []
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            known_bird = 0
	    bird,err = known_birds_f[0],known_birds_f[1]
            for i in range(len(known_birds_f)):
                if (Num.fabs(cand.f-bird)[i] < err[i]):
                    known_bird = 1
                    cand.note = "Freq (%.2f Hz) is within %g Hz " \
                                    "of a known birdie centred at %.2f Hz" % \
                                    (cand.f, err[i], bird[i])
		    bird_removed.append(1)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue
	    bird,err = known_birds_p[0],known_birds_p[1]
            for j in range(len(known_birds_p)):
                if (Num.fabs(cand.p*1000.0-bird)[j] < err[j]):
                    known_bird = 1
                    cand.note = "Period (%.2f ms) is within %g ms " \
                                    "of a known birdie centred at %.2f ms" % \
                                    (cand.f*1000, err[j], bird[j])
		    bird_removed.append(1)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue

    def print_ffacand_summary(self, summaryfilenm=None):
        """Write a summary of all candidates to file (or stdout).

            Input:
                summaryfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if summaryfilenm is None:
            summaryfile = sys.stdout
        elif summaryfilenm in [sys.stdout, sys.stderr]:
            summaryfile = summaryfilenm
        else:
            summaryfile = open(summaryfilenm, "w")
        summaryfile.write("   Candlist contains %d 'good' candidates\n" % \
                            len(self.cands))
        summaryfile.write("      # Known RFI rejects:           %d\n" % \
              len(self.badlists['knownbirds']))
        summaryfile.write("      # Short period rejects:        %d\n" % \
              len(self.badlists['shortperiod']))
        summaryfile.write("      # Long period rejects:         %d\n" % \
              len(self.badlists['longperiod']))
        summaryfile.write("      # Missed threshold:            %d\n" % \
              len(self.badlists['threshold']))
        summaryfile.write("      # Duplicate candidates:        %d\n" % \
              len(self.duplicates))
        summaryfile.write("      # Candidates with DM problems: %d\n" % \
              len(self.badlists['dmproblem']))
        if summaryfilenm not in [None, sys.stdout, sys.stderr]:
            summaryfile.close()
  
    def write_ffacand_report(self, reportfilenm=None):
        """Write a report of all bad candidates to file (or stdout).

            Input:
                reportfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if reportfilenm is None:
            reportfile = sys.stdout
        else:
            reportfile = open(reportfilenm, "w")
        reportfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                       "SNR".center(8) + "P(ms)".center(14) + "numhits".center(9) + "\n")
        badcands = self.get_all_badcands()
        for badcand in badcands:
            reportfile.write("%s (%d)\n" % (str(badcand), len(badcand.hits)))
            reportfile.write("    Note: %s\n\n" % badcand.note)
        if reportfilenm is not None:
            reportfile.close()
	
    def to_file(self, candfilenm=None):
        """Write Candlist to file (or stdout).
            
            Input:
                candfilenm: Name of file to write to. If None,
                    write to stdout. (Default: write to stdout).

            Outputs:
                None
        """
        if candfilenm is None:
            candfile = sys.stdout
        else:
            candfile = open(candfilenm, "w")
        if os.stat(candfilenm).st_size == 0 :
		candfile.write("#" + "file:".center(70)+"candnum".center(40) + \
	"P (ms)".center(40) + "SNR".center(15) + "DM".center(15) + \
                       "dt (ms)".center(15) + 'numhits'.center(15)+"\n")
        for goodcand in self.cands:
            candfile.write("%s \n" % (str(goodcand)+('('+str(len(goodcand.hits))+')').center(15) ))
        if candfilenm is not None:
	    candfile.close()
	    candfile2 = open(candfilenm, "r")
	    next(candfile2)
	    lines = [line.split() for line in candfile2]
	    for l in range(len(lines)):
		lines[l][0] = lines[l][0].split('.dat',1)[0]
		lines[l][0] = lines[l][0] +'.dat:'+str(l+1)+'\t'
		del lines[l][1]
		lines[l][1] = (lines[l][1]).center(40)
		lines[l][2] = (lines[l][2]).center(15)
		lines[l][3] = (lines[l][3]).center(15)
		lines[l][4] = (lines[l][4]).center(15)
		lines[l][5] = (lines[l][5]).center(15)
	    #length = len(lines[-1][0])
	    candfile2 = open(candfilenm, "w")
	    with open(candfilenm, 'w') as fout:
		fout.write("#" + "file:candnum".center(67)+"  P(ms)".center(41) +
		 ("SNR").center(20) +("DM").center(14) + ("dt(ms)").center(18)+ ('numhits').center(12) +'\n')
    		for el in lines:
        		fout.write('{0}\n'.format(' '.join(el)))
	    	
	    fout.close()
	    candfile2.close()
            

#==========================================


def ffa_candlist_from_candfile(filename, trackbad=False, trackdupes=False):
    candfile = open(filename, 'r')
# CHANGED FEW THINGS W.R.T. SIFTING
    # First identify the length of the observation searched
    for line in candfile:
        if line.startswith(" Number of bins in the time series"):
            numsamp = int(line.split()[-1])
        if line.startswith(" Width of each time series bin (sec)"):
            dt = float(line.split()[-1])
	if line.startswith(" Dispersion measure (cm-3 pc)"):
	    dm = float(line.split()[-1])
	    DMstr = str(dm)
    tobs = numsamp * dt
    # Go back to the start of the file to read the candidates
    candfile.seek(0)

    cands = []
    candnums = []
    current_goodcandnum = 0
    final_cands = False
    i = 0
    for line in candfile:
	if line.startswith('#'):
		final_cands = True
		continue
	if (not line.startswith(' ')) and (not line.isspace()):
		split_line = line.split()
	 	candnum = (split_line[0])
        	p = float(split_line[1])       # Spin period in sec
		f = 1.0/p    		       # Spin freq in Hz	
		snr = float(split_line[2])
		dt = float(split_line[3])
		if final_cands :
			filename = str(split_line[0])
			filename = filename.split('.dat',1)[0]
			candnum = i
			p = float(split_line[1])/1000.
			f = 1.0/p
			dm = float(split_line[3])
			dt = float(split_line[4])/1000.
			DMstr = str(dm)
		binn = tobs*f
	    	if DM_re.search(filename) == None and (not final_cands): 
				DMstr = '9999999.99'
            	if (not DM_re.search(filename) == None) and final_cands: 
				DMstr = DM_re.search(filename).groups()[0]
            	cands.append(FFACandidate(candnum,p, snr, dt,binn, dm,DMstr, filename, tobs, final_cands))
		i+=1 
    candfile.close()
    return FFACandlist(cands, trackbad=trackbad, trackdupes=trackdupes)


def ffa_read_candidates(filenms, prelim_reject=True, track=False):
    """Read in accelsearch candidates from the test ACCEL files.
        Return a Candlist object of Candidate instances.

        Inputs:
            filenms: A list of files to read candidates from.
            prelim_reject: If True, perform preliminary rejection of
                candidates. (Default: True)
            track: If True, keep track of bad/duplicate candidates.
                (Default: False)

    """
    candlist = FFACandlist(trackbad=track, trackdupes=track)
    # CHANGED THE READING HERE 
    fi = open(filenms,'r')
    fns=fi.readlines()
    fns=[item.strip() for item in fns]
    numfiles = len(fns)
    ii = 0
    if fns:
        print "\nReading candidates from %d files...." % len(filenms)
        for filenm in fns:
            curr_candlist = ffa_candlist_from_candfile(filenm, trackbad=track, trackdupes=track)
            candlist.extend(curr_candlist)
            sys.stdout.write(" Read %d of %d files (%d cands)\r" % (ii+1, numfiles, len(candlist)))
            sys.stdout.flush()
	    ii+=1
        print "\nDone"
    else:
        print "Error:  There are no candidate files to read!"
    return candlist
