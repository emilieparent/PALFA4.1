import numpy as np
import math
from scipy import *
import scipy.signal
import sys
import re
from operator import itemgetter
import glob
import subprocess 

import ffa_sifting

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

def get_timeseries(beam):
	""" 
	Reads the time series
	Inputs:
		beam : str, file name of the .dat file (dedispersed time series)
	Return :
		ts   : list, time series, dtype = 'float32'
		name : str, name of the beam without the .dat 
	"""
	name = beam
	ts = list(np.fromfile(name,dtype='float32'))
	if name.endswith('.dat'):
		name = re.sub('\.dat$',"",name)
	return ts, name 
	
	
def get_info_beam(beam):
	""" 

	Reads the .inf file, so that the FFA knows what is the duration of the observation, 
	the sampling interval and the DM.
	Inputs:
		beam : str, file name of the .inf file ('.inf' must be included in beam)
	Return :
		T    : float, the lenght of the observation (T in sec)
		dt   : float, the sampling interval (dt in sec)
		DM   : float, dispersion measure	
	"""
	inffile = open(beam,'r')
	for line in inffile:
        	if line.startswith(" Number of bins in the time series"):
           	 	N = int(line.split()[-1])
        	if line.startswith(" Width of each time series bin (sec)"):
            		dt = float(line.split()[-1])
		if line.startswith(" Dispersion measure (cm-3 pc)"):
			DM = float(line.split()[-1])
	T = N * dt
	return T,dt,DM
	
def write_inf(name,file_to_write,write_DM=True):
	"""
	Will write information about the beam at the end of the pre-candidates 
	list (i.e. beam name + "_precands.ffa").
	This is required for sifting (needs the DM).
	Inputs:
		name : str, name of the beam (where it will get the .inf file)
		file_to_write : str, name of the "_precands.ffa" file
		write_DM : bool, write the line telling the dispersion measure or not. 
			   (False when writing final cands)
	"""
	inf_file = open(name + '.inf','r')
	for line in inf_file:
		if (not write_DM) and line.startswith(" Dispersion measure (cm-3 pc)"):
			continue
    		file_to_write.write(line)
	inf_file.close()
	return file_to_write

		
#==================	   Functions related to downsampling 	        ==================


def select_factor(ts, m, mm):
    """ This will delete the last element of ts until it has a factor within [m,mm]
	There is a maximum number of element that you delete in order to get a factor.
	select_factor(ts, m, mm):
	Inputs:
		ts :  array (time series usually)
		m  : minimum factor
		mm : maximum factor
	Returns:
		ts : array (possibly shorter)
		x[0]: first factor that is within the range [minimum,maximum]
    """
    ts = np.array(ts)
    a=np.array(factors(len(ts)))
    x = a[(a >=m) & (a <=mm)]
    counts =0
    while len(x)==0 and counts <50 :
	ts = np.delete(ts,-1)
	a=np.array(factors(len(ts)))
	x = a[(a >=m) & (a <=mm)]
	counts+=1
	if counts >= 50 :
		print "Having a hard time finding a downsampling factor to match the desired time 			resolution. "
		print "Try changing : "+'\n'+'1) the desired time resolution '+'\n'+\
		' 2) by how much you are welling to vary from that time resolution'+'\n'+\
		' 3) How much sample you are willing to delete (default = max 40) '
		sys.exit()

    return ts,x[0] 


def factors(x):
    """	This function takes a number (for FFA purpose, x= lenght of time series) and returns the factors in 		an array, increasing order

	factors(x):
	Inputs :
		x : int
	Returns :
		facts: array of factors of the input x 
    """
    facts = []
    for i in range(1, x + 1):
        if x % i == 0:facts.append(i)
    return facts


def make_factor_10(x):
	"""
	   Deletes the last element of x until x is a factor of 10
	   return array factor of 10. 
	   make_factor_10(x)"
	   Inputs :	   
		x : array
	   Returns:
		x : shortened array (if input is not a factor of 10, otherwise returns input)

	"""
	while len(x)%10!=0:
		x=np.delete(x,-1)
	return x

def normalize(lst):
    """
	normalize using the sum of the points
	normalize(lst):
	Inputs:	
		lst: list, must be 1-D
	Returns 
		Normalized array 
    """
    s = sum(lst)
    return array(map(lambda x: float(x)/s, lst))



def downsample(vector, factor):
    """
       Downsample (i.e. co-add consecutive numbers) a short section
       of a vector by an integer factor. It does not average consecutive bins,
       it adds them.
       downsample(vector, factor)
       Inputs:
		vector: array to be downsampled
		factor: factor by which vector is downsampled (int)
       Returns:
		newvector : downsampled vector
		

    """
    if (len(vector) % factor):
        print "Lenght of 'vector' is not divisible by 'factor'=%d!" % factor
        sys.exit()
    newvector = np.reshape(vector, (len(vector)/factor, factor))
    return np.add.reduce(newvector, 1)


def set_dws(data,downsamp,minimum,maximum):
    """ This will delete last element of the time series so that it can be 
	downsampled by the chosen amount, unless it requires deleting more than 50 elements.
	Inputs :
		data: array (or list), time series
		downsamp: downsampling factor (int)
		minimum:  int, the minimum downsampling factor you are willing to accept
		maximum:  int, the maximum downsampling factor you are willing to accept
    """
    i=0
    data=list(data)
    d = 1
    while d == 1 and i<=60:
    	 data.pop()
   	 n=len(data)
   	 downsamp = select_factor(factors(n),minimum,maximum)
	 i+=1
    if i>60:
	 print'downsamp=1'
	 sys.exit() 												
    return array(data),downsamp	

def forced_dws_2phase(data):
		"""
		Performs downsampling by factor of 2, returns 2 arrays :
		2 different phases (1+2,3+4,.. and 2+3,4+5,..)
		Inputs:
			data: list, time series
		Returns:
			data1, data2 : lists, downsampled time series 
					by factor of 2, 2 phases.
		"""	
		if len(data)%2 == 1:
			data = np.delete(data,-1)
		data1=downsample(data,2)
		f=data[0]
		l=data[-1]
		f_l=(f+l)
		data2=np.delete(data,0)
		data2=np.delete(data2,-1)
		data2=downsample(data2,2)
		data2=np.append(data2,f_l)
		return data1,data2


def forced_dws_3phase(data):
		"""
		Performs downsampling by factor of 3, returns 3 arrays : 
		3 different phases (1+2+3,4+5+6,.. and 2+3+4,5+6+7,.. 
		and 3+4+5,6+7+8,...)
		Inputs:
			data: list, time series
		Returns:
			data1, data2, data3 : lists, downsampled time series 
					       by factor of 3, 3 phases.
		"""	
		if len(data)%3 == 1:
			data = np.delete(data,-1)
		if len(data)%3 == 2:
			data = np.delete(data,-1)
			data = np.delete(data,-1)
		data1=downsample(data,3)
		f=data[0]
		ff=data[1]
		l=data[-1]
		ll=data[-2]
		f_l_ll=(f+l+ll)
		f_ff_l=(f+ff+l)

		data2=np.delete(data,0)
		data2=np.delete(data2,-1)
		data2=np.delete(data2,-1)
		data2=downsample(data2,3)
		data2=np.append(data2,f_l_ll)

		data3=np.delete(data,0)
		data3=np.delete(data3,0)
		data3=np.delete(data3,-1)
		data3=downsample(data3,3)
		data3=np.append(data3,f_ff_l)


		return np.array(data1),np.array(data2),np.array(data2)		


def forced_dws(data,factor):
     """
     Will delete the last elements of a list until the lenght of the 
     list becomes dividible by factor
     Inputs:
	data : list, time series
	factor: int, factor by which to want to downsample the time series 
     Returns:
	data : array, shortened time series (unless already factor of factor)  
     
     """
     while(len(data)%factor!=0):
         data=np.delete(data,-1)
     data=downsample(data,factor)
     return np.array(data)


#==================	   Candidates related functions	        ==================
class ffa_cands(object):
    """
    FFA candidates has 3 parameters: 
    	period (sec)
	SNR  
	dt (sec), the sampling interval  
    """
    def __init__(self):
	self.periods = np.array([])
	self.SNRs =  np.array([])
	self.dts = np.array([])

    def add_cand(self,p,SNR,dt):
	self.periods = np.append(self.periods, p)
	self.SNRs = np.append(self.SNRs, SNR)
	self.dts = np.append(self.dts,dt)

    def print_content(self):
	print "  Period : ", self.periods,  ", SNR: ",\
	self.SNRs, ",dt: ", self.dts
	print "	len(P): ",len(self.periods), ",   len(SNR): ",\
	len(self.SNRs)," ,   len(dts): ",len(self.dts)


cands = ffa_cands()
def cands_to_file(cands,name,suffix):
	"""
	Takes candidates and write it to a file, not the final candidates list.
	Must apply sifting on this list afterward.
	Inputs:
		cands : object that belongs to ffa_cands() class 
		name  : str, name of the beam
	Returns : None
	Writes cands to a file : beam name  + " _precands.ffa" 
	"""
	k=int(1)
	fo = open(name+suffix,'a')
	for i in range(len(cands.periods)):
		fo.write(str(k)+'\t'+'\t'+\
		str(cands.periods[i])+'\t'+'\t'+str(cands.SNRs[i])\
		+'\t'+'\t'+str(cands.dts[i])+'\n')
		k+=1
	print "Wrote ", str(k), "candidates in : "+name+suffix
	write_inf(name,fo)
	fo.close()


def apply_sifting(candsfile, output_name):
	"""
	Apply the sifting on the "_precands.ffa" once done with a ".dat" file.
	It removes duplicate candidates, and removes harmonically related periods.
	Writes the final list of candidates to "_cands.ffa"
	Does not include the final DM sifting ( see siftDM() ).
	Inputs:
		candsfile_list	: a text file that has the name of the "_precands.ffa" to sift
		output_name	: The name you whish to give to the final file (set as the name
				  of the beam + "_cands.ffa" in the ffa script).	
	Returns : None
	Writes the final cands to file
	"""
	candidates = ffa_sifting.ffa_read_candidates(candsfile)
	candidates.remove_duplicate_candidates()
	candidates.remove_harmonics()
	print 'Wrote ',len(candidates),' final candidates in ', output_name
	candidates.to_file(candfilenm = output_name)

	

	

#==================	Signal-to-noise functions	==================

def simple_SNR(folds, sigma, added_profs):
    """ Return a very simple signal-to-noise for a profile.
        Works for narrow duty-cycle since  the S/N is max_value/std
	For each M profiles, returns a value of SNR (i.e, output is a list of lenght M)
    """
    M, P0 = folds.shape
    prof_std = 1.0/(np.sqrt(M-added_profs)*sigma)
    snr = (folds.max(axis=1)-np.median(folds, axis=1))*prof_std
    look_for_nan(snr)
    return snr


def look_for_nan(SNs):
	"""
	For each trial period, makes sure that the S/N is not "nan" 
	Inputs : 
		SNs: list of S/N ratio for this trial period (FFA goes through M different 
		     profiles for each trial period, each M profile has a S/N, and the list
		     SNs is a list of all these M S/N ratios
	Returns : None, 
	Will exit if one of the S/N is nan
	"""
	if SNs.any()=='nan':
		print 'Signal to noise = nan'
		sys.exit()


def param_sn_uniform(SNs):
	"""
	Will make the list of SNs uniform (when considering S/N distributions for 
	different duty cycle tested) 
	param_sn_uniform(SNs):
		SNs : list of signal to noise for all tested period 
	Returns :
		mode : float to be substracted from the list of S/N
		MAD  : float that divides the list of S/N 
	
	"""
	hist = np.histogram(SNs,bins=400,normed=True)
	yval = hist[0]
	xval = np.linspace(hist[1].min(),hist[1].max(),num=400)
	mode = xval[yval.argmax()]
	#sigma = Median absolute deviation, the 1.4826 is a scaling factor for white noise.
	sigma = 1.4826*np.median(np.abs(SNs-np.median(SNs)))
	hist = np.histogram((SNs-mode)/sigma,bins=400,normed=True)
	return  mode, sigma


