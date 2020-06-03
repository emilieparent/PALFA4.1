#!/usr/bin/env python
import numpy as np
import math
import time
import sys
import os
import ast
import re
import matplotlib.pyplot as plt
import scipy.signal
import glob
import subprocess
import optparse

import ffa_tools as ft
import ffa_stages as fs
import FFA_cy as FFA

time_total = time.time()

def func(List,IDs):
    return [List[i] for i in IDs]
    
    
def main():
    """
    Runs the Fast Folding Algorithm on a de-dispersed time series
    The parameters are : the mininum duty-cycle to look for,
                 the Signal-to-Noise treshold to select candidates 
                 the ranges of trial periods (set in config file)
                 the minimum sampling intervals corresponding to 
                         those period ranges (set in config file)
    
    """

    parser = optparse.OptionParser(prog="ffa.py", \
                        version="Emilie Parent (Spring, 2016)", \
                        description="Applies the Fast Folding Algorithm on "\
                " a de-dispersed time series")
    parser.add_option('--mindc',dest='mindc',type = 'float', \
        help="Minimum duty-cycle to look for. Default is 0.5%. "\
            "Options are: 0.5, 1., 1.5. It multiplies the list "\
            "of minimum sampling intervals that has to be tested in each "\
            "subranges of periods by 2 (mindc=1) or 3 (mindc =1.5)",\
            default = None)

    parser.add_option('--SN_tresh',dest='SN_tresh', type='float', \
            help="Signal-to-noise treshold for picking candidates", \
            default = None)

    parser.add_option('--plot',dest='periodograms', \
             action='store_true',  help="Produces periodograms "\
            "for few different duty cycles, default is False")


    options, args = parser.parse_args()
    periodograms = options.periodograms 

    ffa_time=time.time()
    # -------------     Declare list variable    ------------------
    all_SNs_x1,all_SNs1 = [],[]
    all_Ps_x1,all_Ps1 = [] ,[]
    all_SNs_x2,all_SNs_2 = [],[]
    all_Ps_x2,all_Ps_2= [],[]
    all_SNs_x3,all_SNs_3 = [],[]
    all_Ps_x3,all_Ps_3 = [],[]
    SNs1 = []
    SNs2_phase1, SNs2_phase2 = [], []
    SNs4_phase1, SNs4_phase2 ,SNs4_phase3, SNs4_phase4 = [],[],[],[]
    SNs8_phase1 ,SNs8_phase2 ,SNs8_phase3, SNs8_phase4 = [],[],[],[] 
    SNs8_phase5 ,SNs8_phase6 ,SNs8_phase7, SNs8_phase8 = [],[],[],[]
    SNs3_phase1 ,SNs3_phase2, SNs3_phase3 = [], [], []
    SNs9_phase1 ,SNs9_phase2 ,SNs9_phase3 ,SNs9_phase4 = [],[],[],[]
    SNs9_phase5 ,SNs9_phase6 ,SNs9_phase7 ,SNs9_phase8 ,SNs9_phase9 = [],[],[],[],[]
    SNs27_phase1, SNs27_phase2, SNs27_phase3, SNs27_phase4 = [], [], [], []
    SNs27_phase5, SNs27_phase6, SNs27_phase7, SNs27_phase8 = [], [], [], []
    SNs27_phase9, SNs27_phase10, SNs27_phase11, SNs27_phase12 = [], [], [], []
    SNs27_phase13, SNs27_phase14, SNs27_phase15, SNs27_phase16 = [], [], [], []
    SNs27_phase17, SNs27_phase18, SNs27_phase19, SNs27_phase20 = [], [], [], []
    SNs27_phase21, SNs27_phase22, SNs27_phase23, SNs27_phase24 = [], [], [], []
    SNs27_phase25, SNs27_phase26,SNs27_phase27 = [], [], []
    Ps1, Ps2, Ps4, Ps8, Ps3, Ps9, Ps27 = [],[],[],[],[],[],[]
    dts_x1,dts_1 = [],[]
    dts_x2,dts_2 = [],[]
    dts_x3,dts_3 = [],[]
    dt1s,dt2s,dt4s,dt8s,dt3s,dt9s,dt27s = [],[],[],[],[],[],[]
    
    # -------------     Read configuration file        ------------------

    import Config_ffa
    p_ranges = Config_ffa.p_ranges
    dt_list = Config_ffa.dt_list
    if options.SN_tresh ==None: 
        SN_tresh = Config_ffa.SN_tresh
    else: 
        SN_tresh = options.SN_tresh
    if options.mindc ==None: 
        mindc = Config_ffa.mindc
    else:
        mindc = options.mindc
    if periodograms ==None:
         periodograms = False
    
    # -------------     Select beam    ------------------
    beam = sys.argv[1]
    ts,name = ft.get_timeseries(beam)
    T,dt,DM = ft.get_info_beam(name+'.inf')
    N = int(T/dt)
    
    # -------------     Downsampling    ------------------
    # mindc : minimum duty-cycle tested, set in cfg file 
    if mindc > 0.5:
        if mindc ==1:     dt_list = 2*dt_list 
        if mindc ==1.5: dt_list = 3*dt_list
    dwn_ideal = int(dt_list[0]/dt)
    
    # min,max :give a range of accepted downsampling amount - for initial downsampling
    minimum_dwn,maximum_dwn  = dwn_ideal-int(dwn_ideal*0.10),dwn_ideal+int(dwn_ideal*0.20)
    ts,dwn = ft.select_factor(ts,minimum_dwn,maximum_dwn)
    
    print 'Downsampling the time series ..' 
    ts = ft.downsample(ts, dwn)
    
    # -------------     Detrending    ------------------
    
    print "Detrending .."
    #window size : over which statistics are computed 
    window_size = 50*int(len(ts)/T)        
    break_points = np.arange(0,len(ts),window_size) 
    ts = scipy.signal.detrend(ts,bp=break_points)
    
    # Normalize w.r.t. maximum
    ts = ts/max(ts)
    sigma_total= np.std(ts)
    
    #===================================    FFA    ===================================
    print "Entering FFA"
    dt= T/len(ts)
    
    # count_lim: used in stage 2 and 3; how many consecutive downsamplings 
    count_lim = 2        
    
    # Going through subsequent sub-ranges of periods (set in config_ffa.cfg)
    # Each range of periods has it own initial sampling interval
    for num in range(len(p_ranges)):
        if num > 0:
            dwn_ideal = int(dt_list[num]/dt)
            if (dwn_ideal ==1) or (dwn_ideal == 0) :
                dwn_ideal =2 
            minimum_dwn,maximum_dwn = dwn_ideal-int(dwn_ideal*0.05),dwn_ideal+int(dwn_ideal*0.15)
            ts,dwn = ft.select_factor(ts,minimum_dwn,maximum_dwn)
            ts = ft.downsample(ts, dwn)
            sigma_total=sigma_total*np.sqrt(dwn)
            dt = T/len(ts)
        print "  Folding, period range of ", p_ranges[num], " ..."
        all_SNs_x1, all_Ps_x1, dts_x1 = fs.ffa_code_stage1(ts, dt,T, sigma_total,p_ranges[num][0],\
                p_ranges[num][1], count_lim,name)
        all_SNs1.append(all_SNs_x1), all_Ps1.append(all_Ps_x1), dts_1.append(dts_x1)
        all_SNs_x2, all_Ps_x2 , dts_x2 = fs.ffa_code_stage2(ts, dt,T, sigma_total,p_ranges[num][0],\
                p_ranges[num][1], count_lim,name)
        all_SNs_2.append(all_SNs_x2), all_Ps_2.append(all_Ps_x2), dts_2.append(dts_x2)
        all_SNs_x3, all_Ps_x3 , dts_x3 = fs.ffa_code_stage3(ts, dt,T, sigma_total,p_ranges[num][0],\
                p_ranges[num][1], count_lim,name)
        all_SNs_3.append(all_SNs_x3), all_Ps_3.append(all_Ps_x3), dts_3.append(dts_x3)
    
    
    # -------------         end of FFA        ------------------
    
    # Format the lists of S/N, periods, sampling intervals
    for i in range(len(all_Ps1)):
        Ps1.extend(all_Ps1[i]), dt1s.extend(dts_1[i]) , SNs1.extend(all_SNs1[i])
        SNs2_phase1.extend(all_SNs_2[i][0]), SNs2_phase2.extend(all_SNs_2[i][1])
        Ps2.extend(all_Ps_2[i][0]), dt2s.extend(dts_2[i][0])
        SNs4_phase1.extend(all_SNs_2[i][2]), SNs4_phase2.extend(all_SNs_2[i][3])
        SNs4_phase3.extend(all_SNs_2[i][4]), SNs4_phase4.extend(all_SNs_2[i][5])
        Ps4.extend(all_Ps_2[i][1]), dt4s.extend(dts_2[i][1])
        SNs3_phase1.extend(all_SNs_3[i][0]), SNs3_phase2.extend(all_SNs_3[i][1])
        SNs3_phase3.extend(all_SNs_3[i][2])
        Ps3.extend(all_Ps_3[i][0]), dt3s.extend(dts_3[i][0])
        SNs9_phase1.extend(all_SNs_3[i][3]), SNs9_phase2.extend(all_SNs_3[i][4])
        SNs9_phase3.extend(all_SNs_3[i][5]), SNs9_phase4.extend(all_SNs_3[i][6])
        SNs9_phase5.extend(all_SNs_3[i][7]), SNs9_phase6.extend(all_SNs_3[i][8])
        SNs9_phase7.extend(all_SNs_3[i][9]), SNs9_phase8.extend(all_SNs_3[i][10])
        SNs9_phase9.extend(all_SNs_3[i][11])
        Ps9.extend(all_Ps_3[i][1]), dt9s.extend(dts_3[i][1])
        if count_lim ==2:
            SNs8_phase1.extend(all_SNs_2[i][6]), SNs8_phase2.extend(all_SNs_2[i][7])
            SNs8_phase3.extend(all_SNs_2[i][8]), SNs8_phase4.extend(all_SNs_2[i][9])
            SNs8_phase5.extend(all_SNs_2[i][10]), SNs8_phase6.extend(all_SNs_2[i][11])
            SNs8_phase7.extend(all_SNs_2[i][12]), SNs8_phase8.extend(all_SNs_2[i][13])
            Ps8.extend(all_Ps_2[i][2]), dt8s.extend(dts_2[i][2])
            SNs27_phase1.extend(all_SNs_3[i][12]), SNs27_phase2.extend(all_SNs_3[i][13])
            SNs27_phase3.extend(all_SNs_3[i][14]), SNs27_phase4.extend(all_SNs_3[i][15])
            SNs27_phase5.extend(all_SNs_3[i][16]), SNs27_phase6.extend(all_SNs_3[i][17])
            SNs27_phase7.extend(all_SNs_3[i][18]), SNs27_phase8.extend(all_SNs_3[i][19])
            SNs27_phase9.extend(all_SNs_3[i][20]), SNs27_phase10.extend(all_SNs_3[i][21])
            SNs27_phase11.extend(all_SNs_3[i][22]), SNs27_phase12.extend(all_SNs_3[i][23])
            SNs27_phase13.extend(all_SNs_3[i][24]), SNs27_phase14.extend(all_SNs_3[i][25])
            SNs27_phase15.extend(all_SNs_3[i][26]), SNs27_phase16.extend(all_SNs_3[i][27])
            SNs27_phase17.extend(all_SNs_3[i][28]), SNs27_phase18.extend(all_SNs_3[i][29])
            SNs27_phase19.extend(all_SNs_3[i][30]), SNs27_phase20.extend(all_SNs_3[i][31])
            SNs27_phase21.extend(all_SNs_3[i][32]), SNs27_phase22.extend(all_SNs_3[i][33])
            SNs27_phase23.extend(all_SNs_3[i][34]), SNs27_phase24.extend(all_SNs_3[i][35])
            SNs27_phase25.extend(all_SNs_3[i][36]), SNs27_phase26.extend(all_SNs_3[i][37])
            SNs27_phase27.extend(all_SNs_3[i][38])
            Ps27.extend(all_Ps_3[i][2]), dt27s.extend(dts_3[i][2])
    Ps1, SNs1 = np.concatenate(Ps1), np.concatenate(SNs1)

    print "Folding done "
    time_tot =  (time.time() - ffa_time)
    print ( " --- %.7s seconds is the FFA time ---" % time_tot),'\n'


    # ==============================     Post - FFA    =========================================
    # Calculate the mode and the MAD of array of S/N, for 7 (or 5) diff. duty cycles.
    # Assumes that the statistics are the same when duty cycle is the same
    
    loc1, scale1  = ft.param_sn_uniform(SNs1)
    loc2, scale2 = ft.param_sn_uniform(SNs2_phase1)
    loc4, scale4 = ft.param_sn_uniform(SNs4_phase1)
    loc3, scale3 = ft.param_sn_uniform(SNs3_phase1)
    loc9, scale9 = ft.param_sn_uniform(SNs9_phase1)
    if count_lim==2:
        loc8, scale8 = ft.param_sn_uniform(SNs8_phase1)
        loc27, scale27 = ft.param_sn_uniform(SNs27_phase1)
        
    # Making lists of arrays ; all S/Ns for each duty cycle in one object 
    list_SNS = [SNs1,SNs2_phase1, SNs2_phase2, SNs4_phase1, SNs4_phase2, SNs4_phase3, SNs4_phase4, 
            SNs3_phase1, SNs3_phase2, SNs3_phase3, SNs9_phase1, SNs9_phase2, SNs9_phase3,    
            SNs9_phase4, SNs9_phase5, SNs9_phase6, SNs9_phase7, SNs9_phase8, SNs9_phase9]
    
    list_PS = [Ps1, Ps2,Ps2,  Ps4,Ps4,Ps4,Ps4,  Ps3,Ps3,Ps3,  Ps9,Ps9,Ps9,Ps9,Ps9,Ps9,Ps9,Ps9,Ps9]
    list_DTS = [dt1s, dt2s,dt2s, dt4s,dt4s,dt4s,dt4s, dt3s,dt3s,dt3s,   
            dt9s,dt9s,dt9s,dt9s,dt9s,dt9s,dt9s,dt9s,dt9s]
    
    list_locs = [loc1, loc2,loc2, loc4,loc4,loc4,loc4, loc3,loc3,loc3, 
            loc9,loc9,loc9,loc9,loc9,loc9,loc9,loc9,loc9] 
    
    list_scales = [scale1, scale2,scale2, scale4,scale4,scale4,scale4, scale3,scale3,scale3,
            scale9,scale9,scale9,scale9,scale9,scale9,scale9,scale9,scale9]
    
    if count_lim ==2:
        list_SNS = list_SNS + [SNs8_phase1, SNs8_phase2, SNs8_phase3, SNs8_phase4, SNs8_phase5, 
                    SNs8_phase6, SNs8_phase7, SNs8_phase8, SNs27_phase1, SNs27_phase2, 
                    SNs27_phase3, SNs27_phase4, SNs27_phase5, SNs27_phase6, SNs27_phase7, 
                    SNs27_phase8, SNs27_phase9, SNs27_phase10, SNs27_phase11, SNs27_phase12,
                    SNs27_phase13, SNs27_phase14, SNs27_phase15, SNs27_phase16, SNs27_phase17,
                    SNs27_phase18, SNs27_phase19, SNs27_phase20, SNs27_phase21, SNs27_phase22,
                    SNs27_phase23, SNs27_phase24, SNs27_phase25, SNs27_phase26, SNs27_phase27]
        list_PS = list_PS + [Ps8, Ps8, Ps8, Ps8, Ps8, Ps8, Ps8, Ps8,
                    Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, 
                    Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27, Ps27,                     Ps27, Ps27, Ps27 ]
        list_DTS = list_DTS + [dt8s, dt8s, dt8s, dt8s, dt8s, dt8s, dt8s, dt8s, 
                    dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, 
                    dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, 
                    dt27s, dt27s, dt27s, dt27s, dt27s, dt27s, dt27s ]
    
        list_locs = list_locs +  [loc8, loc8, loc8, loc8, loc8, loc8, loc8, loc8, 
                    loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, 
                    loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, loc27, 
                    loc27, loc27, loc27, loc27, loc27, loc27, loc27] 
        list_scales = list_scales + [scale8,  scale8, scale8, scale8, scale8, scale8, scale8, scale8, 
                    scale27, scale27, scale27, scale27, scale27, scale27, scale27, scale27, 
                    scale27, scale27, scale27, scale27, scale27, scale27, scale27, scale27, 
                    scale27, scale27, scale27, scale27, scale27, scale27, scale27, scale27, 
                    scale27, scale27, scale27]
    
    # ======================   Picking Candidates    ==================================
    print "Picking candidates ..."

    # write cands: only True when you get to the end. 
    write_cands, write_cands[-1] = [False]*len(list_SNS), True
    
    # Make all S/Ns uniform
    list_SNS = [(list_SNS[i] -list_locs[i])/list_scales[i] for i in range(len(list_SNS))]
    
    j = [list_SNS[i] >= SN_tresh for i in range(len(list_SNS))]
    good_id = [np.nonzero(j[i]) for i in range(len(j))]
    good_id = np.concatenate(np.array(good_id))
    good_p = [func(list_PS[i],good_id[i]) for i in range(len(good_id))]
    good_sn = [func(list_SNS[i],good_id[i]) for i in range(len(good_id))]
    good_dt = [func(list_DTS[i],good_id[i]) for i in range(len(good_id))]
    good_p  = np.concatenate(good_p)
    good_sn = np.concatenate(good_sn)
    good_dt = np.concatenate(good_dt)
    good_p = [round(good_p[i],4) for i in range(len(good_p))]
    good_sn = [round(good_sn[i],4) for i in range(len(good_sn))]
    good_dt = [round(good_dt[i],5) for i in range(len(good_dt))]

    cands = ft.ffa_cands()
    cands.add_cand(good_p,good_sn,good_dt)
    ft.cands_to_file(cands,name,'_precands.ffa')

    # Check if candidates were selected
    print SN_tresh
    if len(cands.periods)==0:
        print "No cands with S/N > ",SN_tresh," were detected"
        sys.exit()

    print 'Sifting the candidates ...'
    #making precands list (must be sifted)
    candsfile_str = name+'.txt'
    candsfile_list = open(candsfile_str ,'w')
    candsfile_list.write(name+'_precands.ffa')
    candsfile_list.close()
    if os.stat(name+'_precands.ffa').st_size >0:
        ft.apply_sifting(candsfile_str,name+'_cands.ffa')    

    
    print "Completed ", name
    print "Total time for FFA: ",time.time()-time_total        
    subprocess.call(["rm",candsfile_str])
    subprocess.call(["rm",name+'_precands.ffa'])
    
    # Produces periodograms for few duty cycles
    if periodograms :
        plt.figure(figsize=(20,10))
    
        ymin=0
        ymax = np.array([list_SNS[0].max(), list_SNS[1].max(), list_SNS[2].max(), list_SNS[10].max(),\
            list_SNS[12].max(), list_SNS[15].max(), list_SNS[18].max()]).max()
        
        plt.subplot(311)
        plt.suptitle("Periodograms"+'\n'+name)
        # sampling interval : initial
        plt.plot(Ps1,list_SNS[0],color='k',linewidth=1.3,label='No Downsample')
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=0.1,xmax=30)
        plt.ylim(ymin,ymax)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(frameon=False,prop={'size':10})
    
        plt.subplot(312)
        # sampling interval : 2X initial
        plt.plot(Ps2,list_SNS[1],color='blue',linewidth=1.0,label='Downsample x2 (phase 1)')
        plt.plot(Ps2,list_SNS[2],color='steelblue',linewidth=1.0,label='Downsample x2 (phase 2)')
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=0.1,xmax=30)
        plt.ylim(ymin,ymax)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(frameon=False,prop={'size':10})
    
        plt.subplot(313)
        # sampling interval : 9X initial
        plt.plot(Ps9,list_SNS[10],color='r',linewidth=1.0,label='Downsample x9 (phase 1)')
        plt.plot(Ps9,list_SNS[12],color='indianred',linewidth=1.0,label='Downsample x9 (phase 3)')
        plt.plot(Ps9,list_SNS[15],color='tomato',linewidth=1.0,label='Downsample x9 (phase 6)')
        plt.plot(Ps9,list_SNS[18],color='maroon',linewidth=1.0,label='Downsample x9 (phase 9)')
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.xlabel('Period (s)',fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=0.1,xmax=30)
        plt.ylim(ymin,ymax)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(frameon=False,prop={'size':10})

        plt.savefig(name+'.png')
        print "Periodogram available : ", name+'.png'


if __name__=='__main__':
    main()
