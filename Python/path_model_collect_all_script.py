#%matplotlib widget
from multiprocessing import Pool
import sys
import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
sys.path.append( os.path.abspath('/gpfs/users/livermore/jerk_finder/') )
from jerks import jerks
import chaosmagpy as cp

# %%
# Download the spherical harmonic time-series output from the path model
import os.path
if not os.path.exists('Gauss_Bsurf.mat'):
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
# %%
# Run the model on a lat/long grid over the 400 years spanning jerk 9
# Assume 10% error
# Collect all the histograms into a list and save to disk - this takes a while...
# AF attempt at a multithreading approach
jerk_number = 8 # in Python indexing
time_yearly = np.arange(jerk_times[jerk_number]-200,jerk_times[jerk_number]+200+1)

run_components=[0, 1, 2]

ncomp = len(run_components)

SV_error = 10
SV_MIN = -400
SV_MAX = 400
discretise_size = 100
TIMES = time_yearly
TIMES_MIN = TIMES.min()
TIMES_MAX = TIMES.max()
NUM_DATA = len(TIMES)
THIN = 1000
NBINS = 100
credible = 0
RUNNING_MODE = 1
burn_in = 10000


NSAMPLE = 5000 + burn_in

Y_MIN = -400
Y_MAX = 400

K_MIN = 0
K_MAX = 100

sigmas = np.array([20,10,20], dtype = float)
        # sigma_change_value = sigmas(1)
        # sigma_move = sigmas(2)
        # sigma_birth = sigmas(3)
        
# ****************************************
# Define diagnostic time bins here
# ****************************************
time_intervals_edges = np.arange(-210,211,5) + jerk_times[jerk_number]
time_intervals_nbins = time_intervals_edges.shape[0]-1


# ****************************************


CP_hist_save = np.zeros( (len(run_components),time_intervals_nbins), dtype=int )

ntheta = 33
nphi = 72

radius = 6371.2
phis = np.linspace(-180, 180 , nphi, endpoint=False, dtype=float)
thetas = np.linspace(-80,80, ntheta, endpoint=True, dtype=float) + 90.
theta_grid, phi_grid = np.meshgrid(thetas, phis)
thetaphi_g = np.transpose(np.vstack((theta_grid.flatten(), phi_grid.flatten())))
npt = np.shape(thetaphi_g)[0]

def my_calc_par(thetaphi_g):

    npt = np.shape(thetaphi_g)[0]
    loc_results = []
    theta_l = thetaphi_g[0]
    phi_l = thetaphi_g[1]
    Br, Btheta, Bphi = cp.model_utils.synth_values(coeffs, radius, theta_l, phi_l, nmax=13)
    Br_yearly,Btheta_yearly, Bphi_yearly = np.interp(time_yearly, time, Br ), np.interp(time_yearly, time, Btheta ), np.interp(time_yearly, time, Bphi )
    Bx_dot, By_dot, Bz_dot = -np.gradient(Btheta_yearly,time_yearly), np.gradient(Bphi_yearly,time_yearly), -np.gradient(Br_yearly,time_yearly)
    for i in run_components:
        #print(i)
        if i == 0:
            SV = Bx_dot
        elif i == 1:
            SV = By_dot
        else:
            SV = Bz_dot

        delta_SV = 0.01 * SV_error * (SV.max() - SV.min()) * np.ones(NUM_DATA,dtype=float)
        
        TIME_grid = np.linspace(TIMES_MIN, TIMES_MAX, discretise_size, endpoint=True)
      

        
        Acceptance_rates=np.zeros(4)
        AV = np.zeros(discretise_size,dtype=float)
        SUP = np.zeros(discretise_size,dtype=float)
        INF = np.zeros(discretise_size,dtype=float)
        MEDIAN = np.zeros(discretise_size,dtype=float)
        MODE = np.zeros(discretise_size,dtype=float)
        CP_hist_run = np.zeros( time_intervals_nbins, dtype=int )
        MARGINAL_DENSITY = np.zeros( (discretise_size,NBINS),dtype=float )
        N_CP_hist = np.zeros( K_MAX, dtype=int)
        av_delta = np.zeros( time_intervals_nbins,dtype=float )
        median_delta = np.zeros( time_intervals_nbins,dtype=float )
        upper_delta = np.zeros( time_intervals_nbins,dtype=float )
        lower_delta = np.zeros( time_intervals_nbins,dtype=float )

        (Acceptance_rates, SUP, INF, AV, MEDIAN, MODE, CP, av_delta,median_delta, upper_delta, lower_delta,MARGINAL_DENSITY, N_CP_hist) = jerks.rjmcmc(
         sigmas=sigmas, burn_in=burn_in,
         nsample=NSAMPLE, num_data=NUM_DATA, times=TIMES, y=SV, delta_y=delta_SV,
         y_min=SV_MIN, y_max=SV_MAX, times_min=TIMES_MIN, times_max=TIMES_MAX, k_min=K_MIN,
         k_max=K_MAX, discretise_size=discretise_size, time_intervals_nbins = time_intervals_nbins,time_intervals_edges = time_intervals_edges,
         thin=THIN, nbins=NBINS, credible=credible, running_mode=RUNNING_MODE)

        # save the model
        fac = (NSAMPLE-burn_in)/THIN
        loc_results += [theta_l, phi_l,i,CP_hist_run[:]/fac,AV,av_delta,median_delta, Acceptance_rates]
    return loc_results


nproc = int(os.getenv("OMP_NUM_THREADS"))

if __name__ == '__main__':
    if nproc > 1:
        results = []
        with Pool(nproc) as p:
            loc_results = p.map(my_calc_par, thetaphi_g)
            results += loc_results
    else:
        results = []
        for ipt in range(npt):
            loc_results = my_calc_par( thetaphi_g[ipt])
            results += [loc_results]

import pickle
with open("Jerk{0:d}_5x5_{1:d}M_size{2:d}_error{3:d}_5yr.results".format(jerk_number+1,int((NSAMPLE-burn_in)/1e6),discretise_size,SV_error), "wb") as fp:   #Pickling
    pickle.dump([TIMES_MAX, TIMES_MIN, K_MIN, K_MAX, THIN, burn_in,sigmas,NSAMPLE, discretise_size, NUM_DATA,NBINS, RUNNING_MODE],fp)
    pickle.dump(results, fp)

