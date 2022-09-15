
import sys
import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
sys.path.append( os.path.abspath('..') )
from jerks import jerks
import chaosmagpy as cp

# %%
# Download the spherical harmonic time-series output from the path model
import os.path
if not os.path.exists('./Gauss_Bsurf.mat'):
    os.system('wget --no-verbose http://morpho.ipgp.fr/4DEarth/Gauss_Bsurf.mat')

# %%
# import the dataset
import h5py
filepath = 'Gauss_Bsurf.mat'
arrays = {}
f = h5py.File(filepath,'r')
for k, v in f.items():
    arrays[k] = np.array(v)

# %%
coeffs = arrays['gnm'][:,:].T
time = arrays['timers'].flatten()
print( 'Shape of gmn array: ', arrays['gnm'].shape )
print ('Times (in years) of output', time )

# %%
# jerk times as defined in catalogue
jerk_times = [4600,5750,2920, 1915, 6490,7300,7620,7840,8880,9673,10590,12620,13411,13546]
"""
# %%
"""


radius = 6371.0
theta = 80
phi = 30
jerk_number = 2 # in Python indexing

# sample every year
time_yearly = np.arange(jerk_times[jerk_number]-200,jerk_times[jerk_number]+200+1)

Br, Btheta, Bphi = cp.model_utils.synth_values(coeffs, radius, theta, phi,nmax=13)
Br_yearly,Btheta_yearly, Bphi_yearly = np.interp(time_yearly, time, Br ), np.interp(time_yearly, time, Btheta ), np.interp(time_yearly, time, Bphi )
Bx_dot, By_dot, Bz_dot = -np.gradient(Btheta_yearly,time_yearly), np.gradient(Bphi_yearly,time_yearly), -np.gradient(Br_yearly,time_yearly)


SV_error = 10

discretise_size = 100
TIMES = time_yearly
NUM_DATA = len(TIMES)
TIMES_MIN = TIMES.min()
TIMES_MAX = TIMES.max()
CP_NBINS = int(TIMES_MAX - TIMES_MIN) #one per year
    
time_intervals_edges = np.arange(-210,211,20) + jerk_times[jerk_number]
time_intervals_nbins = time_intervals_edges.shape[0]-1


    
   
SV = By_dot


delta_SV = SV_error * 0.01 * (SV.max() - SV.min()) * np.ones(NUM_DATA,dtype=float)
SV_MIN = -400
SV_MAX = 400

K_MIN = 0
K_MAX = 100
sigmas = np.array([10,5,10],dtype = float)
TIME_grid = np.linspace(TIMES_MIN, TIMES_MAX, discretise_size, endpoint=True)
# sigma_change_value = sigmas(1)
# sigma_move = sigmas(2)
# sigma_birth = sigmas(3)

THIN = 1000
NBINS = 100
credible = 0.0
RUNNING_MODE = 1
burn_in = 10000
Acceptance_rates=np.zeros(4)
AV = np.zeros(discretise_size,dtype=float)
SUP = np.zeros(discretise_size,dtype=float)
INF = np.zeros(discretise_size,dtype=float)
MEDIAN = np.zeros(discretise_size,dtype=float)
MODE = np.zeros(discretise_size,dtype=float)
CP_hist_run = np.zeros( CP_NBINS, dtype=int )
MARGINAL_DENSITY = np.zeros( (discretise_size,NBINS),dtype=float )
N_CP_hist = np.zeros( K_MAX, dtype=int)
delta_slope = np.zeros( time_intervals_nbins,dtype=float )

print( theta ,phi )
NSAMPLE = 2000000+burn_in
(Acceptance_rates, SUP, INF, AV, MEDIAN, MODE, CP_hist_run, delta_slope, MARGINAL_DENSITY, N_CP_hist) = jerks.rjmcmc(
 sigmas=sigmas, burn_in=burn_in,
     nsample=NSAMPLE, num_data=NUM_DATA, times=TIMES, y=SV, delta_y=delta_SV,
     y_min=SV_MIN, y_max=SV_MAX, times_min=TIMES_MIN, times_max=TIMES_MAX, k_min=K_MIN,
     k_max=K_MAX, discretise_size=discretise_size, time_intervals_nbins = time_intervals_nbins,time_intervals_edges = time_intervals_edges,
     thin=THIN, nbins=NBINS, credible=credible, running_mode=RUNNING_MODE)
CP_hist_run = CP_hist_run / ( (NSAMPLE-burn_in)/THIN )
print(radius, NSAMPLE-burn_in)
print(delta_slope)
print( Acceptance_rates )
