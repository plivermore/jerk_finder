# %%
#%matplotlib widget
from multiprocessing import Pool
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
import chaosmagpy as cp

# adjust the path according to the computer on which the code is running:
import socket
print('Computer: ' + socket.gethostname())

if socket.gethostname()[0:4] == 'ncpu':
    sys.path.append( os.path.abspath('/gpfs/users/livermore/jerk_finder/') )
else:
    sys.path.append( os.path.abspath('/Users/earpwl/Library/CloudStorage/OneDrive-UniversityofLeeds/home/codes/Github/jerk_finder') ) #This line needs to point to whereever the jerks module has been compiled.

from jerks import jerks
from find_jerks import find_jerks

TEST_RUN = False


# import the dataset


model = cp.load_CHAOS_matfile('/Users/earpwl/Library/CloudStorage/OneDrive-UniversityofLeeds/home/models/CHAOS-7/CHAOS-7.18.mat')  # load model


time = np.linspace(1999, 2024, 3*(2024-1999)+1)  # decimal years
mjd = cp.data_utils.dyear_to_mjd(time, leap_year=False)

nproc = int(os.getenv("OMP_NUM_THREADS"))


# Assume 3nT/yr error
SV_error = 3.0

# components:
run_components=[0,1,2]


radius = 6371.2

# Run the model on a lat/long grid
phis = np.arange(-180,180,5,dtype=float)
thetas = np.arange(10,175,5,dtype=float)

if TEST_RUN:
    phis = np.arange(-180,180,20,dtype=float)
    thetas = np.arange(10,175,20,dtype=float)
    print('Running as test run...')
theta_grid, phi_grid = np.meshgrid(thetas, phis)
thetaphi_g = np.transpose(np.vstack((theta_grid.flatten(), phi_grid.flatten())))
npt = np.shape(thetaphi_g)[0]



       
K_MIN = 0
K_MAX = 100
THIN = 100
NBINS = 100
credible = 0
RUNNING_MODE = 1
burn_in = 10000
NSAMPLE = 1000000+burn_in

if TEST_RUN:
    burn_in = 100
    NSAMPLE = 1000+burn_in
    
sigmas = np.array([3.0,1.0,3.0])

def my_calc_par(thetaphi_g):

    npt = np.shape(thetaphi_g)[0]
    loc_results = []
    theta_l = thetaphi_g[0]
    phi_l = thetaphi_g[1]
    
    for cmpt in run_components:
    
        dBr, dBt, dBp = model.synth_values_tdep(mjd, radius, theta_l, phi_l, nmax=13, deriv=1)

        if cmpt == 0:
            SV = -dBt
        elif cmpt == 1:
            SV = dBp
        else:
            SV = -dBr
            
        SV_MIN = SV.min()
        SV_MAX = SV.max()

        discretise_size = 100
        NUM_DATA = len(time)
        TIMES_MIN = time.min()
        TIMES_MAX = time.max()
        delta_SV = SV_error * np.ones(NUM_DATA,dtype=float)

        
        # ****************************************
        
        size_jerk_data = 0
        jerk_data = np.zeros( ( K_MAX *(NSAMPLE-burn_in)//THIN,2),dtype=float )

        TIME_grid = np.linspace(TIMES_MIN, TIMES_MAX, discretise_size, endpoint=True)

        Acceptance_rates=np.zeros(4)
        AV = np.zeros(discretise_size,dtype=float)
        SUP = np.zeros(discretise_size,dtype=float)
        INF = np.zeros(discretise_size,dtype=float)
        MEDIAN = np.zeros(discretise_size,dtype=float)
        MODE = np.zeros(discretise_size,dtype=float)
        MARGINAL_DENSITY = np.zeros( (discretise_size,NBINS),dtype=float )
        N_CP_hist = np.zeros( K_MAX+1, dtype=int)  #indices 0...KMAX
            
        #saved history for animations
        jerk_history = np.zeros( ((NSAMPLE-burn_in)//THIN,2), dtype=int)
        model_history = np.zeros( (discretise_size, (NSAMPLE-burn_in)//THIN),dtype=float)

        (Acceptance_rates, SUP, INF, AV, MEDIAN, MODE,
        MARGINAL_DENSITY, N_CP_hist, jerk_data, size_jerk_data, \
        model_history, jerk_history) = jerks.rjmcmc(
         sigmas=sigmas, burn_in=burn_in,
         nsample=NSAMPLE, num_data=NUM_DATA, times=time, y=SV, delta_y=delta_SV,
         y_min=SV_MIN, y_max=SV_MAX, times_min=TIMES_MIN, times_max=TIMES_MAX, k_min=K_MIN,
         k_max=K_MAX, discretise_size=discretise_size,
         thin=THIN, nbins=NBINS, credible=credible, running_mode=RUNNING_MODE)


        data_timing = jerk_data[:,0]
        data_amplitude = jerk_data[:,1]
        
        
        number_time_bins = NUM_DATA - 2 #ignore edge data
        time_range = [TIMES_MIN+1, TIMES_MAX-1]  #exclude end effects
        jerks_info = find_jerks(data_timing, data_amplitude, number_time_bins, time_range)

    
        loc_results += [theta_l, phi_l,cmpt,jerks_info]
    return loc_results


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
            results += loc_results


    import pickle
    with open("Jerks_CHAOS718_abserror_3_5x5_{0:1d}M.results".format((NSAMPLE-burn_in)//1000000), "wb") as fp:   #Pickling
    # export meta data
    
        meta = [SV_error, K_MIN, K_MAX, THIN, burn_in, NSAMPLE, sigmas, thetas, phis]
    
        pickle.dump(meta, fp)
        pickle.dump(results, fp)
