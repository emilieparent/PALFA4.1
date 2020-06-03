#!/usr/bin/env python
import sys
import subprocess 
import numpy as np
import re
import glob
import optparse

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import ffa_tools as ft
import ffa_stages as fs
import ffa_sifting 

"""
Run this script once you have applied the FFA on all de-dispersed time series.
The candidate file inputed should include candidates from multiple DMs (if it applies)

You need to give it a  '_rfifind.inf' file you want to write at the end of the final candidate list
"""


def get_zaplist(zapfn):
    fctr, width = np.loadtxt(zapfn, usecols=(-2,-1), unpack=True)
    return fctr, width


def final_sifting_ffa(basenm, candfile_list, output_file, zapfn=[]):
	fout = open(output_file,'w')
	fout.write("#" + "file".center(20)+"		                P(ms)".center(40) +
		 "SNR".center(24) +"DM".center(4) + "dt (ms)".center(44) +'\n')

	#take all the candidates from the list of files in candfile_list
	#writes them in output_file
	for fil in candfile_list:
		temp = open(fil,'r')
		lines_temp = [line.split() for line in temp]
		del lines_temp[0]
		for l in lines_temp:
			string = l[0].center(20)+"	"+l[1].center(40) +l[2].center(2) +\
					l[3].center(24) + l[4].center(24) +'\n'
			fout.write(string)
	#writes information at the end of output_file
	ft.write_inf(basenm+'_rfifind',fout,write_DM=False)
	fout.close()
	# sifting need a file that contains the list of files to sift 
	# only one file in this file, in our case
	# candfile_for_sifting has the file on which sifting must be applied 

	candfile_for_sifting = 'for_sifting_'+output_file
	ftemp = open(candfile_for_sifting,'w')
	ftemp.write(output_file)
	ftemp.close()

	#sifting
	dmlist= [re.sub('\_cands.ffa$',"",candfile_list[i]).split('_DM',1)[-1] \
				      for i in range(len(candfile_list))]
	while any(".dat" in dmlist[i] for i in range(len(dmlist))):
		dmlist=[re.sub('\.dat$','',dmlist[i]) for i in range(len(dmlist))]
	candidates = ffa_sifting.ffa_read_candidates(candfile_for_sifting)
	if zapfn != None:
		known_birds_f = tuple(get_zaplist(zapfn))
		candidates.reject_knownbirds(known_birds_f)
	candidates.remove_duplicate_candidates()
	candidates.remove_DM_problems(numdms=2, dmlist=dmlist, low_DM_cutoff = 2.0 )
	try: 
		candidates.remove_harmonics()
	except: 
		print "No more candidates to sift."
    #### make summary plot for FFA cands ####
        candidates.plot_ffa_summary(usefreqs=False)	
        plt.title("%s FFA Summary" % basenm)
        plt.savefig(basenm+".ffacands.summary.png")	
    #### make rejects plot for FFA cands ####
        candidates.plot_ffa_rejects(usefreqs=False)	
        plt.title("%s FFA Rejected Cands" % basenm)
        plt.savefig(basenm+".ffacands.rejects.png")	
    #### write reports and print summary ####
        candidates.print_ffacand_summary(basenm+".ffacands.summary")    
        candidates.write_ffacand_report(basenm+".ffacands.report")
	candidates.to_file(candfilenm = output_file)

	#candfile_for_sifting is useless after sifting; delete it
	subprocess.call(["rm",candfile_for_sifting])
	finf = open(output_file, 'a')
	ft.write_inf(basenm+'_rfifind',finf,write_DM=False)
	print 'Wrote ',len(candidates),' final cands in ', output_file
	return candidates


def main():
	parser = optparse.OptionParser(prog="ffa_final.py", \
                        version="Emilie Parent (Spring, 2016)", \
                        description="Performs the final sifting of FFA "\
				                        "candidates (at all DMs), and writes "\
				                        "a .txt file containing the final candidates.")
	parser.add_option('--zaplist',dest='zaplist',type = 'string', \
		                    help="Zaplist file produced by zaplist (.txt file). "\
		                          "Contains frequencies to be zapped while "\
		                          "sifting the FFA candidates")


	options, args = parser.parse_args()
	basenm_inf = sys.argv[1]
	basenm = re.sub('\_rfifind.inf$',"",basenm_inf)
	candfile_list = glob.glob(basenm+"*_DM*_cands.ffa")
	zapfn = options.zaplist
	output_file = basenm+'_cands.ffa'
	ffa_cands = final_sifting_ffa(basenm,candfile_list,output_file, zapfn)
	
	
	
if __name__=='__main__':
    main()

