# %%
#%matplotlib widget
from multiprocessing import Pool
import sys
import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
#sys.path.append( os.path.abspath('..') )
sys.path.append( os.path.abspath('/gpfs/users/livermore/jerk_finder/') )
from jerks import jerks
import chaosmagpy as cp

# Find jerks:

# 1. Find peaks in marginal histogram of changepoints satisfying criteria like min distance between jerks and jerk prominence,
# 2. Find the subset of peaks for which the associated modal amplitude (from the 2d histogram) exceeds a certain threshold, label the peak a a jerk
# 3. Find the jerk width by seeing how far either side of the peak the 1d histogram drops to tau of the value at the peak.

# defaults:
#
def find_jerks(data_timing, data_amplitude, number_time_bins, time_range,
    number_amplitude_bins = 100, amplitude_range = [-50,50],
    jerk_min_amplitude = 1, min_distance = 4, \
    peak_prominence = 1, tau = 0.14):
    '''
    Function: find_jerks
    
    Returns a list of jerk attributes
    
        Parameters:
            data_timing (float, array): changepoint times
            data_amplitude: (float, array): changepoint slope changes
            number_time_bins (integer): number of bins in time
            time_range (list): start and end times
            
            Optional:
            jerk_min_amplitude (float): threshold for jerk identification
            min_distance (float): minimum distance between jerks
            prominence (float): prominence of jerks
            tau (float): threshold to characterise width of jerks
            
        Returns a list of jerk attributes in the form
            [jerk0, jerk1,...]
            where (for example) jerk0 is [jerk_time, jerk_amplitude, min_jerk_time, max_jerk_time]
            All four values are float; jerk_time is the most likely time for a jerk of amplitude jerk_amplitude, but it is uncertain with bounds given by min_jerk_time and max_jerk_time.
            
    Algorithm:
    1. Find peaks in marginal histogram of changepoints satisfying: min distance between jerks and jerk prominence.
    2. Find the subset of peaks for which the associated slope change of maximum probability (from the 2d histogram) exceeds the threshold of jerk_min_amplitude.  See Scipy documentation for definition of prominence.
    3. For these peaks ("jerks"), the jerk width is quantified by identifying the smallest time window that defines a drop in the 1d histogram of changepoint timing by tau on both sides.
    If the changepoint histogram is Gaussian in shape, then 2 sigma corresponds to a drop in height of e^-2, or about 0.14. This defines then (roughly) a 95% credible interval.
    

    '''
    
    jerks_overall_info = []

    

    counts, xedges, yedges = np.histogram2d(data_timing, data_amplitude, bins=(number_time_bins,number_amplitude_bins), range=[time_range, amplitude_range])
    marginal_counts, marginal_xedges = np.histogram( data_timing, bins=number_time_bins, range=time_range)
    from scipy import signal

    peaks = signal.find_peaks(marginal_counts,distance=min_distance, prominence = peak_prominence * marginal_counts.mean() )

    for peak_index in peaks[0]:
        jerk_info = np.zeros(4)
        index = counts[peak_index,:].argmax()  #find index corresponding to max probability
        jerk_amplitude = 0.5 * ( yedges[index] + yedges[index+1] )  #use the centre point of the bin
        if abs(jerk_amplitude) > jerk_min_amplitude:
            jerk_info[0] = 0.5 * ( marginal_xedges[peak_index] + marginal_xedges[peak_index+1] )  #use the centre point of the bin
            jerk_info[1] = jerk_amplitude

            min_jerk_time = marginal_xedges[0]
            max_jerk_time = marginal_xedges[-1]

            for i in range(peak_index,marginal_counts.shape[0]):  #count upwards in time to find when the jerk window ends
                if marginal_counts[i] < marginal_counts[peak_index] * tau:
                    max_jerk_time = marginal_xedges[i]  #i here indicates the first bin which lies outside the jerk. The edge of the jerk is thus taken to be the left edge of this bin.
                    break

            for i in range(peak_index,0,-1):  # count downwards to find when the jerk window ends
                if marginal_counts[i] < marginal_counts[peak_index] * tau:
                    min_jerk_time = marginal_xedges[i+1]  #i here indicates the first bin which lies outside the jerk. The edge of the jerk is thus taken to be the right edge of this bin.
                    break

            jerk_info[2], jerk_info[3] = min_jerk_time, max_jerk_time

            jerks_overall_info.append(jerk_info)
    return jerks_overall_info

    

save_all_data = True


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

# define window lengths
window_length = 400
overlap = 50

#coeffs = coeffs[time <= 8000, :]
#time = time[time <= 8000]


nproc = int(os.getenv("OMP_NUM_THREADS"))

window_start = np.arange(int(time.min()), int(time.max()), window_length)
window_end= window_start + window_length

window_overlap_start = [max(i-overlap,int(time.min()))  for i in window_start]
window_overlap_end = [min(i+overlap,int(time.max()))  for i in window_end]

print('Number of windows defined: ', len(window_start))
print('Window list:', list(zip(window_start, window_end)))


    
# %%
# Run the model on a lat/long grid
# Assume 10% error

# components:
run_components=[0,1,2]




radius = 6371.2

phis = np.arange(-180,180,5,dtype=float)
thetas = np.arange(10,175,5,dtype=float)


theta_grid, phi_grid = np.meshgrid(thetas, phis)
thetaphi_g = np.transpose(np.vstack((theta_grid.flatten(), phi_grid.flatten())))
npt = np.shape(thetaphi_g)[0]

#print('Theta grid: ', thetas)
#print('Phi grid: ', phis)
#print('Overall: {0:d} grid points '.format(len(theta_grid.flatten() )))

# relative errors
SV_error = 10
       
K_MIN = 0
K_MAX = 100
THIN = 100
NBINS = 100
credible = 0
RUNNING_MODE = 1
burn_in = 10000
NSAMPLE = 1000000+burn_in

#burn_in = 100
#NSAMPLE = 1000+burn_in
relative_sigmas = np.array([0.08, 0.02, 0.08])
        
def my_calc_par(thetaphi_g):

    npt = np.shape(thetaphi_g)[0]
    loc_results = []
    theta_l = thetaphi_g[0]
    phi_l = thetaphi_g[1]
    
    size_jerk_data_stitch  = np.zeros( len(run_components) )
    jerk_data_stitch0, jerk_data_stitch1, jerk_data_stitch2 = np.empty(shape=[0, 2]), np.empty(shape=[0, 2]), np.empty(shape=[0, 2])
    
    


    for cmpt in run_components:
    
        # Find min/max SV for whole timeseries
        Br, Btheta, Bphi = cp.model_utils.synth_values(coeffs, radius, theta_l, phi_l,nmax=13)
        time_yearly = np.arange(int(time.min()), int(time.max())+1)
        Br_yearly,Btheta_yearly, Bphi_yearly = np.interp(time_yearly, time, Br ), \
        np.interp(time_yearly, time, Btheta ), np.interp(time_yearly, time, Bphi )
        Bx_dot, By_dot, Bz_dot = -np.gradient(Btheta_yearly,time_yearly), \
        np.gradient(Bphi_yearly,time_yearly), -np.gradient(Br_yearly,time_yearly)

        if cmpt == 0:
            SV = Bx_dot
        elif cmpt == 1:
            SV = By_dot
        else:
            SV = Bz_dot
            
        SV_MIN = SV.min()
        SV_MAX = SV.max()

        for window_index in range(len(window_overlap_start)):
        
        # resample SV within time window
            time_yearly = np.arange(window_overlap_start[window_index], window_overlap_end[window_index]+1)
            Br_yearly,Btheta_yearly, Bphi_yearly = np.interp(time_yearly, time, Br ), np.interp(time_yearly, time, Btheta ), np.interp(time_yearly, time, Bphi )
            Bx_dot, By_dot, Bz_dot = -np.gradient(Btheta_yearly,time_yearly), np.gradient(Bphi_yearly,time_yearly), -np.gradient(Br_yearly,time_yearly)
        
            if cmpt == 0:
                SV = Bx_dot
            elif cmpt == 1:
                SV = By_dot
            else:
                SV = Bz_dot

            discretise_size = 100
            TIMES = time_yearly
            NUM_DATA = len(TIMES)
            TIMES_MIN = TIMES.min()
            TIMES_MAX = TIMES.max()
            delta_SV = SV_error * 0.01 * (SV_MAX - SV_MIN) * np.ones(NUM_DATA,dtype=float)

            sigmas = np.array([ (SV_MAX - SV_MIN)*relative_sigmas[0],\
                               (TIMES_MAX - TIMES_MIN)*relative_sigmas[1],\
                               (SV_MAX - SV_MIN)*relative_sigmas[2]],dtype = float)

       
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
            

            (Acceptance_rates, SUP, INF, AV, MEDIAN, MODE,
            MARGINAL_DENSITY, N_CP_hist, jerk_data, size_jerk_data) = jerks.rjmcmc(
             sigmas=sigmas, burn_in=burn_in,
             nsample=NSAMPLE, num_data=NUM_DATA, times=TIMES, y=SV, delta_y=delta_SV,
             y_min=SV_MIN, y_max=SV_MAX, times_min=TIMES_MIN, times_max=TIMES_MAX, k_min=K_MIN,
             k_max=K_MAX, discretise_size=discretise_size,
             thin=THIN, nbins=NBINS, credible=credible, running_mode=RUNNING_MODE)

# save only the part within the actual (non-overlapping) time window:
            mask = [i for i in range(jerk_data.shape[0]) if (i < size_jerk_data) and \
                (jerk_data[i,0] < window_end[window_index]) and \
                (jerk_data[i,0] > window_start[window_index])]
# don't worry about jerks falling exactly on the window boundaries, as this happens with probability 0.

            if cmpt == 0:
                jerk_data_stitch0 = np.vstack((jerk_data_stitch0[:,:],jerk_data[mask,:]))
            elif cmpt == 1:
                jerk_data_stitch1 = np.vstack((jerk_data_stitch1[:,:],jerk_data[mask,:]))
            else:
                jerk_data_stitch2 = np.vstack((jerk_data_stitch2[:,:],jerk_data[mask,:]))
                
            size_jerk_data_stitch[cmpt] += len(mask)

            #print('Theta,phi {5:3.1f} {6:3.1f}; window index: {0:d}; acceptance rates: {1:3.1f}% {2:3.1f}% {3:3.1f}% #{4:3.1f}%'.\
                #  format(window_index, Acceptance_rates[0],Acceptance_rates[1],Acceptance_rates[2],Acceptance_rates[3],\
                #  theta_l, phi_l ))


        # save the model after all the windowed calculations
        

        if cmpt == 0:
            data_timing, data_amplitude = jerk_data_stitch0[:,0], jerk_data_stitch0[:,1]
        elif cmpt==1:
            data_timing, data_amplitude = jerk_data_stitch1[:,0], jerk_data_stitch1[:,1]
        else:
            data_timing, data_amplitude = jerk_data_stitch2[:,0], jerk_data_stitch2[:,1]
            
        TIMES_MAX, TIMES_MIN = int( time.max() ), int( time.min() )
        number_time_bins = (TIMES_MAX - TIMES_MIN) * 2 #one bin per 6 months.
        number_amplitude_bins = 100
        time_range = [TIMES_MIN, TIMES_MAX]
        amplitude_range = [-50,50]
        jerks_info = find_jerks(data_timing, data_amplitude, number_time_bins, time_range, tau = 0.60)

        if save_all_data:
            import pickle
            path = 'Jerk_data_files'
            if not os.path.exists(path):
                os.makedirs(path)
            with open(os.path.join(path,"Theta={0:3.1f}phi={1:3.1f}_{2:1d}M_data.results".format(theta_l, phi_l,(NSAMPLE-burn_in)//1000000)), "wb") as fp:   #Pickling
                pickle.dump([jerk_data_stitch0,jerk_data_stitch1, jerk_data_stitch2], fp)
    
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
    with open("Jerks_5x5_{0:1d}M.results".format((NSAMPLE-burn_in)//1000000), "wb") as fp:   #Pickling
    # export meta data
    
        meta = [SV_error, K_MIN, K_MAX, THIN, burn_in, NSAMPLE, \
        relative_sigmas, window_start, window_end, window_length, \
        overlap, thetas, phis]
    
        pickle.dump(meta, fp)
        pickle.dump(results, fp)
