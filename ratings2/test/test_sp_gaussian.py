#!/usr/bin/python

import sys
import pprint

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import Greys

import sp_utils

import candidate
import utils
from sp_rating_classes import gaussian

def main():
    spdfn = sys.argv[1]

    spd = sp_utils.spd(spdfn)

    cand = candidate.SPCandidate(spd.best_dm, \
                        spd.ra_deg, spd.dec_deg, \
                        spdfn)

    print "Loaded %s" % cand.spdfn
    print "-"*20
    print "    Best DM (cm^-3/pc): %f" % cand.dm
    print "    RA (J2000 - deg): %f" % cand.raj_deg
    print "    Dec (J2000 - deg): %f" % cand.decj_deg

    gauss.add_data(cand)

    print "-"*20
    fit_and_plot(cand, spd)
 
def fit_and_plot(cand, spd):
    data = cand.profile
    n = len(data)
    rms = np.std(data[(n/2):])
    xs = np.linspace(0.0, 1.0, n, endpoint=False)
    G = gauss._compute_data(cand)
    print "    Reduced chi-squared: %f" % (G.get_chisqr(data) / G.get_dof(n))
    print "    Baseline rms: %f" % rms
    print "    %s" % G.components[0]

    fig1 = plt.figure(figsize=(10,10))
    plt.subplots_adjust(wspace=0, hspace=0)

    # upper
    ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
    ax1.plot(xs, data/rms, color="black", label="data")
    ax1.plot(xs, G.components[0].make_gaussian(n), color="red", label="best fit")

    # lower
    ax2 = plt.subplot2grid((3,1), (2,0), sharex=ax1)
    ax2.plot(xs, data/rms - G.components[0].make_gaussian(n), color="black", label="residuals")
    ax2.set_xlabel("Fraction of pulse window")

    plt.figure()
    plt.pcolormesh(xs, spd.waterfall_freq_axis(), spd.data_zerodm_dedisp, cmap=Greys)
    plt.xlabel("Fraction of pulse window")
    plt.ylabel("Frequency (MHz)")
    plt.xlim(0, 1)
    plt.ylim(spd.min_freq, spd.max_freq)

    plt.show()

if __name__=='__main__':
    gauss = gaussian.GaussianProfileClass()
    main()
