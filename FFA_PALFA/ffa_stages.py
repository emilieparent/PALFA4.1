import ffa_tools as f
import FFA_cy as FFA
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from scipy import *
import scipy.signal
import inspect
import sys 
import os
"""
Make histogram
check for why sigma isnt sig*sqrt{dwn}
"""
def ffa_code_stage1(data ,dt , T , sigma_total,p_min, p_max, count_lim, name):
	"""
	ffa_code_stage1 (data , dt , T, period , N , p_min , p_max , count_lim , name):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)

	Returns 
	"""
	# --------------------	   FFA Stage 1	----------------------
	fill_value = np.median(data)
	N = T/dt
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	SNs1 = []
	all_Ps1 = []
	for p0 in P0s:
		p0 = int(p0)
		M_real = float(float(N)/p0)
		added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
		if p0==0 or p0 ==1:
			print '                  It tried to fold with period = 0 or 1 bin'
			continue
		p0_sec = p0*dt
		xwrap = FFA.XWrap2(data,p0, fill_value=fill_value, pow2=True)
		folds = FFA.FFA(xwrap)
		M = folds.shape[0]
		SN = f.simple_SNR(folds, np.std(data),added_profs)  
		P = p0 + (np.arange(M, dtype=np.float) / (M-1))
		Psec=P*dt
		SNs1.append(SN)
		all_Ps1.extend(Psec)
	SNs1 = np.concatenate(SNs1)
	dts = [dt]*len(SNs1)
	return np.array([SNs1]),np.array([all_Ps1]), dts

#______________________________________________________________________


# --------------------	   FFA Stage 2	------------------------------
# -------------------- Extra downsampling : 2 	----------------------

def ffa_code_stage2(data ,dt , T , sigma_total,p_min, p_max, count_lim, name):
	"""
	ffa_code_stage2 (data , dt , T, period , N , p_min , p_max , count_lim , name):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)

	Returns 

	"""	
	fill_value = np.median(data)
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	count=0
	while count<=count_lim:
		if count==0:	
			data_1,data_2 = f.forced_dws_2phase(data)
			new_dt2 = dt*2
			dt2 = new_dt2
			all_Ps2 = []
			N = T/(new_dt2)
			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			time_stage2_1=time.time()
			SNs2_1, SNs2_2 = [], []
			for p0 in P0s2:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin','\n'
					continue
				p0_sec = p0*new_dt2
				time_stage2_1=time.time()
				xwrap_1 = FFA.XWrap2(data_1,p0, fill_value=fill_value, pow2=True)
				xwrap_2 = FFA.XWrap2(data_2,p0, fill_value=fill_value, pow2=True)
				folds_1 = FFA.FFA(xwrap_1)
				folds_2 = FFA.FFA(xwrap_2)	
				SN_1 = f.simple_SNR(folds_1, np.std(data_1), added_profs)
				SN_2 = f.simple_SNR(folds_2, np.std(data_2) , added_profs)

				M = folds_1.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2 

				SNs2_1.append(SN_1)
				SNs2_2.append(SN_2)
				all_Ps2.extend(Psec)
			SNs2_1=np.concatenate(SNs2_1)
			SNs2_2=np.concatenate(SNs2_2)
			dt2s = [dt2]*len(SNs2_1)
				
		if count == 1:
			data_11, data_12 = f.forced_dws_2phase(data_1)
			data_21, data_22 = f.forced_dws_2phase(data_2)	
			new_dt2 = new_dt2*2
			dt4 = new_dt2
			N = T/(new_dt2)
			all_Ps4 = []
			SNs4_11, SNs4_12, SNs4_21, SNs4_22 = [],[],[],[]
			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			for p0 in P0s2:
				p0=int(p0)
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin'
					continue
				p0_sec = p0*new_dt2
				xwrap_11 = FFA.XWrap2(data_11,p0, fill_value=fill_value, pow2=True)
				xwrap_12 = FFA.XWrap2(data_12,p0, fill_value=fill_value, pow2=True)
				xwrap_21 = FFA.XWrap2(data_21,p0, fill_value=fill_value, pow2=True)
				xwrap_22 = FFA.XWrap2(data_22,p0, fill_value=fill_value, pow2=True)
		
				folds_11 = FFA.FFA(xwrap_11)
				folds_12 = FFA.FFA(xwrap_12)
				folds_21 = FFA.FFA(xwrap_21)
				folds_22 = FFA.FFA(xwrap_22)

				SN_11 = f.simple_SNR(folds_11, np.std(data_11) , added_profs)
				SN_12 = f.simple_SNR(folds_12, np.std(data_12) , added_profs)
				SN_21 = f.simple_SNR(folds_21, np.std(data_21) , added_profs)
				SN_22 = f.simple_SNR(folds_22, np.std(data_22) , added_profs)

				M = folds_11.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2 
				SNs4_11.append(SN_11)
				SNs4_12.append(SN_12)
				SNs4_21.append(SN_21)
				SNs4_22.append(SN_22)

				all_Ps4.extend(Psec)
			SNs4_11=np.concatenate(SNs4_11)
			SNs4_12=np.concatenate(SNs4_12)
			SNs4_21=np.concatenate(SNs4_21)
			SNs4_22=np.concatenate(SNs4_22)
			dt4s = [dt4]*len(SNs4_11)

		if count == 2:

			data_111, data_112 = f.forced_dws_2phase(data_11)
			data_121, data_122 = f.forced_dws_2phase(data_12)
			data_211, data_212 = f.forced_dws_2phase(data_21)
			data_221, data_222 = f.forced_dws_2phase(data_22)	
			new_dt2 = new_dt2*2
			dt8 = new_dt2
			all_Ps8 = []
			N = T/(new_dt2)
			SNs8_111 , SNs8_112 ,SNs8_121 ,SNs8_122 =[],[],[],[] 
			SNs8_211 ,SNs8_212 ,SNs8_221 ,SNs8_222 =[],[],[],[]
			P0_start, P0_end = np.floor(p_min/new_dt2), np.ceil(p_max/new_dt2)
			P0s2=np.arange(P0_start,P0_end,1)
			for p0 in P0s2:
				p0=int(p0)

				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real

				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin'
					continue
				p0_sec = p0*new_dt2

				xwrap_111 = FFA.XWrap2(data_111,p0, fill_value=fill_value, pow2=True)
				xwrap_121 = FFA.XWrap2(data_121,p0, fill_value=fill_value, pow2=True)
				xwrap_211 = FFA.XWrap2(data_211,p0, fill_value=fill_value, pow2=True)
				xwrap_221 = FFA.XWrap2(data_221,p0, fill_value=fill_value, pow2=True)
				xwrap_112 = FFA.XWrap2(data_112,p0, fill_value=fill_value, pow2=True)
				xwrap_122 = FFA.XWrap2(data_122,p0, fill_value=fill_value, pow2=True)
				xwrap_212 = FFA.XWrap2(data_212,p0, fill_value=fill_value, pow2=True)
				xwrap_222 = FFA.XWrap2(data_222,p0, fill_value=fill_value, pow2=True)
		
				folds_111 = FFA.FFA(xwrap_111)
				folds_121 = FFA.FFA(xwrap_121)
				folds_211 = FFA.FFA(xwrap_211)
				folds_221 = FFA.FFA(xwrap_221)
				folds_112 = FFA.FFA(xwrap_112)
				folds_122 = FFA.FFA(xwrap_122)
				folds_212 = FFA.FFA(xwrap_212)
				folds_222 = FFA.FFA(xwrap_222)


				SN_111 = f.simple_SNR(folds_111, np.std(data_111) , added_profs)
				SN_121 = f.simple_SNR(folds_121, np.std(data_121) , added_profs)
				SN_211 = f.simple_SNR(folds_211, np.std(data_211) , added_profs)
				SN_221 = f.simple_SNR(folds_221, np.std(data_221) , added_profs)
				SN_112 = f.simple_SNR(folds_112, np.std(data_112) , added_profs)
				SN_122 = f.simple_SNR(folds_122, np.std(data_122) , added_profs)
				SN_212 = f.simple_SNR(folds_212, np.std(data_212) , added_profs)
				SN_222 = f.simple_SNR(folds_222, np.std(data_222) , added_profs)				
				M = folds_111.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt2

				SNs8_111.append(SN_111)
				SNs8_121.append(SN_121)
				SNs8_211.append(SN_211) 
				SNs8_221.append(SN_221) 
				SNs8_112.append(SN_112) 
				SNs8_122.append(SN_122) 
				SNs8_212.append(SN_212) 
				SNs8_222.append(SN_222)	

				all_Ps8.extend(Psec)

			SNs8_111 = np.concatenate(SNs8_111)
			SNs8_121 = np.concatenate(SNs8_121)
			SNs8_211 = np.concatenate(SNs8_211)
			SNs8_221 = np.concatenate(SNs8_221)
			SNs8_112 = np.concatenate(SNs8_112)
			SNs8_122 = np.concatenate(SNs8_122)
			SNs8_212 = np.concatenate(SNs8_212)
			SNs8_222 = np.concatenate(SNs8_222)
			dt8s = [dt8]*len(SNs8_111)

		count+=1
	if count_lim == 1:
		return np.array([SNs2_1,SNs2_2, SNs4_11 , SNs4_12 , SNs4_21 , SNs4_22]), np.array([all_Ps2, all_Ps4]), [dt2s,dt4s]
	if count_lim == 2:
		return np.array([SNs2_1,SNs2_2, SNs4_11 , SNs4_12 , SNs4_21 , SNs4_22, SNs8_111 , SNs8_121 , SNs8_211 , SNs8_221 , SNs8_112 , SNs8_122 , SNs8_212 , SNs8_222 ]), np.array([all_Ps2, all_Ps4, all_Ps8]), [dt2s,dt4s,dt8s]

	
	
# --------------------	   FFA Stage 3	----------------------
# -------------------- -                Extra downsampling : 3 	----------------------
def ffa_code_stage3(data ,dt ,T,sigma_total,p_min,p_max, count_lim,name):
	"""
	ffa_code_stage3 (data , dt , T, period , N , p_min , p_max , count_lim , name):
		- data		:  Time series 
		- dt		:  Sampling interval (s)
		- T		:  Total observative time (s)
		- p_min		:  Minimum period in the subset of trial periods (s)
		- p_max		:  Maximum period in the subset of trial periods (s)
		- count_lim	:  int 1 or 2, used in stage2 and stage 3
				   if count_lim =1, goes to 4*dt and 9*dt 
				   if count_lim =2, goes to 8*dt and 27*dt
		- name		:  Name of the beam (without the extension)

	Returns 
	"""	
	fill_value = np.median(data)
	P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
	P0s = np.arange(P0_start,P0_end,1)
	time_ffa3=time.time()
	count=0
	plot_num=311
	larges_dt3=[]
	w=int(3)
	while count<=count_lim:
		if count==0:	
			data_1,data_2,data_3=f.forced_dws_3phase(data)
			new_dt3 = dt*3
			dt3 = new_dt3
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			all_Ps3 = []
			SNs3_1, SNs3_2, SNs3_3 = [], [], []
			for p0 in P0s3:
				p0=int(p0)
				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin'
					continue
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				p0_sec = p0*new_dt3
				xwrap_1 = FFA.XWrap2(data_1,p0, fill_value=fill_value, pow2=True)
				xwrap_2 = FFA.XWrap2(data_2,p0, fill_value=fill_value, pow2=True)
				xwrap_3 = FFA.XWrap2(data_3,p0, fill_value=fill_value, pow2=True)
			
				folds_1 = FFA.FFA(xwrap_1)
				folds_2 = FFA.FFA(xwrap_2)
				folds_3 = FFA.FFA(xwrap_3)

				SN_1 = f.simple_SNR(folds_1, np.std(data_1) , added_profs)
				SN_2 = f.simple_SNR(folds_2, np.std(data_2) , added_profs)
				SN_3 = f.simple_SNR(folds_3, np.std(data_3) , added_profs)

				M = folds_1.shape[0]

				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 
				
				all_Ps3.extend(Psec)

				SNs3_1.append(SN_1)
				SNs3_2.append(SN_2)
				SNs3_3.append(SN_3)
			
			SNs3_1 = np.concatenate(SNs3_1)
			SNs3_2 = np.concatenate(SNs3_2)
			SNs3_3 = np.concatenate(SNs3_3)
			dt3s = [dt3]*len(SNs3_1)
			
		if count==1:
			w=int(w*3)
			new_dt3 = new_dt3*3
			dt9 = new_dt3	
			all_Ps9 = []
			data_11, data_12, data_13 = f.forced_dws_3phase(data_1)
			data_21, data_22, data_23 = f.forced_dws_3phase(data_2)
			data_31, data_32, data_33 = f.forced_dws_3phase(data_3)
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			SNs9_11 ,SNs9_12 ,SNs9_13 ,SNs9_21 ,SNs9_22 = [],[],[],[],[] 
			SNs9_23 ,SNs9_31 ,SNs9_32 ,SNs9_33  = [],[],[],[]
			for p0 in P0s3:
				p0=int(p0)
				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin'
					continue
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				p0_sec = p0*new_dt3
				xwrap_11 = FFA.XWrap2(data_11,p0, fill_value=fill_value, pow2=True)
				xwrap_21 = FFA.XWrap2(data_21,p0, fill_value=fill_value, pow2=True)
				xwrap_31 = FFA.XWrap2(data_31,p0, fill_value=fill_value, pow2=True)
				xwrap_12 = FFA.XWrap2(data_12,p0, fill_value=fill_value, pow2=True)
				xwrap_22 = FFA.XWrap2(data_22,p0, fill_value=fill_value, pow2=True)
				xwrap_32 = FFA.XWrap2(data_32,p0, fill_value=fill_value, pow2=True)
				xwrap_13 = FFA.XWrap2(data_13,p0, fill_value=fill_value, pow2=True)
				xwrap_23 = FFA.XWrap2(data_23,p0, fill_value=fill_value, pow2=True)
				xwrap_33 = FFA.XWrap2(data_33,p0, fill_value=fill_value, pow2=True)
			
				folds_11 = FFA.FFA(xwrap_11)
				folds_21 = FFA.FFA(xwrap_21)
				folds_31 = FFA.FFA(xwrap_31)
				folds_12 = FFA.FFA(xwrap_12)
				folds_22 = FFA.FFA(xwrap_22)
				folds_32 = FFA.FFA(xwrap_32)
				folds_13 = FFA.FFA(xwrap_13)
				folds_23 = FFA.FFA(xwrap_23)
				folds_33 = FFA.FFA(xwrap_33)

				SN_11 = f.simple_SNR(folds_11, np.std(data_11) , added_profs)
				SN_21 = f.simple_SNR(folds_21, np.std(data_21) , added_profs)
				SN_31 = f.simple_SNR(folds_31, np.std(data_31) , added_profs)
				SN_12 = f.simple_SNR(folds_12, np.std(data_12) , added_profs)
				SN_22 = f.simple_SNR(folds_22, np.std(data_22) , added_profs)
				SN_32 = f.simple_SNR(folds_32, np.std(data_32) , added_profs)
				SN_13 = f.simple_SNR(folds_13, np.std(data_13) , added_profs)
				SN_23 = f.simple_SNR(folds_23, np.std(data_23) , added_profs)
				SN_33 = f.simple_SNR(folds_33, np.std(data_33) , added_profs)


				M = folds_11.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 
				SNs9_11.append(SN_11)
				SNs9_12.append(SN_12)
				SNs9_13.append(SN_13)
				SNs9_21.append(SN_21)
				SNs9_22.append(SN_22)
				SNs9_23.append(SN_23)
				SNs9_31.append(SN_31)
				SNs9_32.append(SN_32)
				SNs9_33.append(SN_33)

				all_Ps9.extend(Psec)
			try:
				SNs9_11 = np.concatenate(SNs9_11)
				SNs9_12 = np.concatenate(SNs9_12)
				SNs9_13 = np.concatenate(SNs9_13)
				SNs9_21 = np.concatenate(SNs9_21)
				SNs9_22 = np.concatenate(SNs9_22)
				SNs9_23 = np.concatenate(SNs9_23)
				SNs9_31 = np.concatenate(SNs9_31)
				SNs9_32 = np.concatenate(SNs9_32)
				SNs9_33 = np.concatenate(SNs9_33)
			except :
				continue	
			dt9s = [dt9]*len(SNs9_11)

		if count==2 :
			new_dt3 = new_dt3*3
			dt27 = new_dt3
			w=int(w*3)	
			data_111, data_112, data_113 = f.forced_dws_3phase(data_11)
			data_211, data_212, data_213 = f.forced_dws_3phase(data_21)
			data_311, data_312, data_313 = f.forced_dws_3phase(data_31)
			data_121, data_122, data_123 = f.forced_dws_3phase(data_12)
			data_221, data_222, data_223 = f.forced_dws_3phase(data_22)
			data_321, data_322, data_323 = f.forced_dws_3phase(data_32)
			data_131, data_132, data_133 = f.forced_dws_3phase(data_13)
			data_231, data_232, data_233 = f.forced_dws_3phase(data_23)
			data_331, data_332, data_333 = f.forced_dws_3phase(data_33)
			
			all_Ps27 = []
			N = T/(new_dt3)
			P0_start, P0_end = np.floor(p_min/new_dt3), np.ceil(p_max/new_dt3)
			P0s3=np.arange(P0_start,P0_end,1)
			SNs27_111, SNs27_211, SNs27_311, SNs27_121 = [], [], [], []
			SNs27_221, SNs27_321, SNs27_131, SNs27_231 = [], [], [], []
			SNs27_331, SNs27_112, SNs27_212, SNs27_312 = [], [], [], []
			SNs27_122, SNs27_222, SNs27_322, SNs27_132 = [], [], [], []
			SNs27_232, SNs27_332, SNs27_113, SNs27_213 = [], [], [], []
			SNs27_313, SNs27_123, SNs27_223, SNs27_323 = [], [], [], []
			SNs27_133,SNs27_233,SNs27_333 = [], [], []

			for p0 in P0s3:
				p0=int(p0)
				if p0==0 or p0 ==1:
					print '                  It tried to fold with period = 0 or 1 bin'
					continue
				M_real = float(float(N)/p0)
				added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
				p0_sec = p0*new_dt3
				xwrap_111 = FFA.XWrap2(data_111,p0, fill_value=fill_value, pow2=True)
				xwrap_112 = FFA.XWrap2(data_112,p0, fill_value=fill_value, pow2=True)
				xwrap_113 = FFA.XWrap2(data_113,p0, fill_value=fill_value, pow2=True)	
				xwrap_121 = FFA.XWrap2(data_121,p0, fill_value=fill_value, pow2=True)
				xwrap_122 = FFA.XWrap2(data_122,p0, fill_value=fill_value, pow2=True)
				xwrap_123 = FFA.XWrap2(data_123,p0, fill_value=fill_value, pow2=True)
				xwrap_131 = FFA.XWrap2(data_131,p0, fill_value=fill_value, pow2=True)
				xwrap_132 = FFA.XWrap2(data_132,p0, fill_value=fill_value, pow2=True)
				xwrap_133 = FFA.XWrap2(data_133,p0, fill_value=fill_value, pow2=True)
				xwrap_211 = FFA.XWrap2(data_211,p0, fill_value=fill_value, pow2=True)
				xwrap_212 = FFA.XWrap2(data_212,p0, fill_value=fill_value, pow2=True)
				xwrap_213 = FFA.XWrap2(data_213,p0, fill_value=fill_value, pow2=True)
				xwrap_221 = FFA.XWrap2(data_221,p0, fill_value=fill_value, pow2=True)
				xwrap_222 = FFA.XWrap2(data_222,p0, fill_value=fill_value, pow2=True)
				xwrap_223 = FFA.XWrap2(data_223,p0, fill_value=fill_value, pow2=True)
				xwrap_231 = FFA.XWrap2(data_231,p0, fill_value=fill_value, pow2=True)
				xwrap_232 = FFA.XWrap2(data_232,p0, fill_value=fill_value, pow2=True)
				xwrap_233 = FFA.XWrap2(data_233,p0, fill_value=fill_value, pow2=True)
				xwrap_311 = FFA.XWrap2(data_311,p0, fill_value=fill_value, pow2=True)
				xwrap_312 = FFA.XWrap2(data_312,p0, fill_value=fill_value, pow2=True)
				xwrap_313 = FFA.XWrap2(data_313,p0, fill_value=fill_value, pow2=True)
				xwrap_321 = FFA.XWrap2(data_321,p0, fill_value=fill_value, pow2=True)
				xwrap_322 = FFA.XWrap2(data_322,p0, fill_value=fill_value, pow2=True)
				xwrap_323 = FFA.XWrap2(data_323,p0, fill_value=fill_value, pow2=True)
				xwrap_331 = FFA.XWrap2(data_331,p0, fill_value=fill_value, pow2=True)
				xwrap_332 = FFA.XWrap2(data_332,p0, fill_value=fill_value, pow2=True)
				xwrap_333 = FFA.XWrap2(data_333,p0, fill_value=fill_value, pow2=True)
				
				folds_111 = FFA.FFA(xwrap_111)
				folds_112 = FFA.FFA(xwrap_112)
				folds_113 = FFA.FFA(xwrap_113)
				folds_211 = FFA.FFA(xwrap_211)
				folds_212 = FFA.FFA(xwrap_212)
				folds_213 = FFA.FFA(xwrap_213)
				folds_311 = FFA.FFA(xwrap_311)
				folds_312 = FFA.FFA(xwrap_312)
				folds_313 = FFA.FFA(xwrap_313)
				folds_121 = FFA.FFA(xwrap_121)
				folds_122 = FFA.FFA(xwrap_122)
				folds_123 = FFA.FFA(xwrap_123)
				folds_221 = FFA.FFA(xwrap_221)
				folds_222 = FFA.FFA(xwrap_222)
				folds_223 = FFA.FFA(xwrap_223)
				folds_321 = FFA.FFA(xwrap_321)
				folds_322 = FFA.FFA(xwrap_322)
				folds_323 = FFA.FFA(xwrap_323)
				folds_131 = FFA.FFA(xwrap_131)
				folds_132 = FFA.FFA(xwrap_132)
				folds_133 = FFA.FFA(xwrap_133)
				folds_231 = FFA.FFA(xwrap_231)
				folds_232 = FFA.FFA(xwrap_232)
				folds_233 = FFA.FFA(xwrap_233)
				folds_331 = FFA.FFA(xwrap_331)
				folds_332 = FFA.FFA(xwrap_332)
				folds_333 = FFA.FFA(xwrap_333)
			

				SN_111 = f.simple_SNR(folds_111, np.std(data_111) , added_profs)
				SN_211 = f.simple_SNR(folds_211, np.std(data_211) , added_profs)
				SN_311 = f.simple_SNR(folds_311, np.std(data_311) , added_profs)
				SN_121 = f.simple_SNR(folds_121, np.std(data_121) , added_profs)
				SN_221 = f.simple_SNR(folds_221, np.std(data_221) , added_profs)
				SN_321 = f.simple_SNR(folds_321, np.std(data_321) , added_profs)
				SN_131 = f.simple_SNR(folds_131, np.std(data_131) , added_profs)
				SN_231 = f.simple_SNR(folds_231, np.std(data_231) , added_profs)
				SN_331 = f.simple_SNR(folds_331, np.std(data_331) , added_profs)
				SN_112 = f.simple_SNR(folds_112, np.std(data_112) , added_profs)
				SN_212 = f.simple_SNR(folds_212, np.std(data_212) , added_profs)
				SN_312 = f.simple_SNR(folds_312, np.std(data_312) , added_profs)
				SN_122 = f.simple_SNR(folds_122, np.std(data_122) , added_profs)
				SN_222 = f.simple_SNR(folds_222, np.std(data_222) , added_profs)
				SN_322 = f.simple_SNR(folds_322, np.std(data_322) , added_profs)
				SN_132 = f.simple_SNR(folds_132, np.std(data_132) , added_profs)
				SN_232 = f.simple_SNR(folds_232, np.std(data_232) , added_profs)
				SN_332 = f.simple_SNR(folds_332, np.std(data_332) , added_profs)
				SN_113 = f.simple_SNR(folds_113, np.std(data_113) , added_profs)
				SN_213 = f.simple_SNR(folds_213, np.std(data_213) , added_profs)
				SN_313 = f.simple_SNR(folds_313, np.std(data_313) , added_profs)
				SN_123 = f.simple_SNR(folds_123, np.std(data_123) , added_profs)
				SN_223 = f.simple_SNR(folds_223, np.std(data_223) , added_profs)
				SN_323 = f.simple_SNR(folds_323, np.std(data_323) , added_profs)
				SN_133 = f.simple_SNR(folds_133, np.std(data_133) , added_profs)
				SN_233 = f.simple_SNR(folds_233, np.std(data_233) , added_profs)
				SN_333 = f.simple_SNR(folds_333, np.std(data_333) , added_profs)
				M = folds_111.shape[0]
				P = p0 + (np.arange(M, dtype=np.float) / (M-1))
				Psec = P *new_dt3 

				SNs27_111.append(SN_111)
				SNs27_211.append(SN_211)
				SNs27_311.append(SN_311)
				SNs27_121.append(SN_121)
				SNs27_221.append(SN_221)
				SNs27_321.append(SN_321)
				SNs27_131.append(SN_131)
				SNs27_231.append(SN_231)
				SNs27_331.append(SN_331)
				SNs27_112.append(SN_112)
				SNs27_212.append(SN_212)
				SNs27_312.append(SN_312)
				SNs27_122.append(SN_122)
				SNs27_222.append(SN_222)
				SNs27_322.append(SN_322)
				SNs27_132.append(SN_132)
				SNs27_232.append(SN_232)
				SNs27_332.append(SN_332)
				SNs27_113.append(SN_113) 
				SNs27_213.append(SN_213)
				SNs27_313.append(SN_313)
				SNs27_123.append(SN_123) 
				SNs27_223.append(SN_223) 
				SNs27_323.append(SN_323) 
				SNs27_133.append(SN_133) 
				SNs27_233.append(SN_233) 
				SNs27_333.append(SN_333) 

				all_Ps27.extend(Psec)

			SNs27_111 = np.concatenate(SNs27_111)
			SNs27_211 = np.concatenate(SNs27_211)
			SNs27_311 = np.concatenate(SNs27_311)
			SNs27_121 = np.concatenate(SNs27_121)
			SNs27_221 = np.concatenate(SNs27_221)
			SNs27_321 = np.concatenate(SNs27_321)
			SNs27_131 = np.concatenate(SNs27_131)
			SNs27_231 = np.concatenate(SNs27_231)
			SNs27_331 = np.concatenate(SNs27_331)
			SNs27_112 = np.concatenate(SNs27_112)
			SNs27_212 = np.concatenate(SNs27_212)
			SNs27_312 = np.concatenate(SNs27_312)
			SNs27_122 = np.concatenate(SNs27_122)
			SNs27_222 = np.concatenate(SNs27_222)
			SNs27_322 = np.concatenate(SNs27_322)
			SNs27_132 = np.concatenate(SNs27_132)
			SNs27_232 = np.concatenate(SNs27_232)
			SNs27_332 = np.concatenate(SNs27_332)
			SNs27_113 = np.concatenate(SNs27_113)
			SNs27_213 = np.concatenate(SNs27_213)
			SNs27_313 = np.concatenate(SNs27_313)
			SNs27_123 = np.concatenate(SNs27_123)
			SNs27_223 = np.concatenate(SNs27_223)
			SNs27_323 = np.concatenate(SNs27_323)
			SNs27_133 = np.concatenate(SNs27_133)
			SNs27_233 = np.concatenate(SNs27_233)
			SNs27_333 = np.concatenate(SNs27_333)
			dt27s = [dt27]*len(SNs27_111)

		count+=1
	if count_lim ==1:
		return np.array([SNs3_1 , SNs3_2 , SNs3_3 , SNs9_11, SNs9_12, SNs9_13, SNs9_21, SNs9_22, SNs9_23, SNs9_31, SNs9_32,SNs9_33]), np.array([all_Ps3,all_Ps9]), [dt3s,dt9s]
	if count_lim ==2:
 		return np.array([SNs3_1 , SNs3_2 , SNs3_3 , SNs9_11, SNs9_12, SNs9_13, SNs9_21, SNs9_22, SNs9_23, SNs9_31, SNs9_32,SNs9_33, SNs27_111 ,SNs27_211 ,SNs27_311 ,SNs27_121 ,SNs27_221 ,SNs27_321 ,
SNs27_131 ,SNs27_231 ,SNs27_331 ,SNs27_112 ,SNs27_212 ,SNs27_312 ,SNs27_122 ,SNs27_222 ,
SNs27_322 ,SNs27_132 ,SNs27_232 ,SNs27_332 ,SNs27_113 ,SNs27_213 ,SNs27_313 ,SNs27_123 ,
SNs27_223 ,SNs27_323 ,SNs27_133 ,SNs27_233 ,SNs27_333 ]), np.array([all_Ps3,all_Ps9, all_Ps27]), [dt3s,dt9s,dt27s]

