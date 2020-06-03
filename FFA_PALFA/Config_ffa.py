
import numpy as np

#configuration file FFA

#[FFA_settings]
p_ranges = np.array([[0.5,1.],[1.,2.],[2.,5.],[5.,10.],[10.,15.],[15,30]])
dt_list = np.array([0.005,0.01,0.02,0.050,0.075,0.120])
SN_tresh = 5

# min_dc = 0.5, 1 or 1.5 
mindc = 0.5
numdms = 2
