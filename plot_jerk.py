import numpy as np
import pickle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

jerk_number = 8 
# jerk times as defined in catalogue
jerk_times = [4600,5750,2920, 1915, 6490,7300,7620,7840,8880,9673,10590,12620,13411,13546]
# sample every year
time_yearly = np.arange(jerk_times[jerk_number]-200,jerk_times[jerk_number]+200+1)
t0 = jerk_times[jerk_number]

filename = "Jerk"+str(jerk_number+1)+"_5x5_20M.results"
with open(filename, "rb") as fp:   # Unpickling
     results = pickle.load(fp)

# make plot of number of jerks over t0 +/- 30 years

height_threshold = 0.3
distance_threshold = 3

x_theta, x_phi,x_count = [],[],[]
y_theta, y_phi,y_count = [],[],[]
z_theta, z_phi,z_count = [],[],[]

for j in range(len(results)):
    theta = results[j][0]
    phi = results[j][1]
    component = results[j][2]
    CP = results[j][3]
    if component == 0:
        x_theta.append(theta); x_phi.append(phi)
        peaks,_ = find_peaks( CP, height = height_threshold, distance = distance_threshold)
        x_count.append(len( [time_yearly[i] for i in peaks if (time_yearly[i] < t0+30 and time_yearly[i] > t0-30)]))

    if component == 1:
        y_theta.append(theta); y_phi.append(phi)
        peaks,_ = find_peaks( CP, height = height_threshold, distance = distance_threshold)
        y_count.append(len( [time_yearly[i] for i in peaks if (time_yearly[i] < t0+30 and time_yearly[i] > t0-30)]))
    
    if component == 2:
        z_theta.append(theta); z_phi.append(phi)
        peaks,_ = find_peaks( CP, height = height_threshold, distance = distance_threshold)
        z_count.append(len( [time_yearly[i] for i in peaks if (time_yearly[i] < t0+30 and time_yearly[i] > t0-30)]))
    t = len( [time_yearly[i] for i in peaks 
                 if (time_yearly[i] < t0+30 and time_yearly[i] > t0-30)])
    if t >= 3:
                 print(str(t) + ' jerks found at theta = ' + str(theta) + ' phi = ' + str(phi) + ' component = ' + str(component))
max_count = 4#max(max(x_count),max(y_count),max(z_count))
# force this to be 4.

cmap = plt.get_cmap('rainbow', max_count-0+1)
print( 'Max jerk count is ' + str(max_count))   
plt.figure()
axes = [0,0,0]
f, (axes[0],axes[1],axes[2]) = plt.subplots(nrows=3, ncols=1, figsize=(7,12),sharex=True, subplot_kw={'projection': ccrs.PlateCarree() }) 

marker_size = 10. 

for i in range(3):
    if i == 0:
        cax = axes[i].scatter(x_phi,90.-np.array(x_theta), s = marker_size, c=x_count,cmap=cmap, vmin=0-0.5, vmax=(max_count)+0.5)
        axes[i].set_title(r'$dB_X/dt, \qquad \Sigma_{{total}} = {0:d},\; \Sigma_{{nz}} = {1:d}$'.format(sum(x_count),np.sum(np.array(x_count)>0)))
        gl = axes[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
    elif i == 1:
        cax = axes[i].scatter(y_phi,90.-np.array(y_theta), s = marker_size, c=y_count,cmap=cmap, vmin=0-0.5, vmax=(max_count)+0.5)
        axes[i].set_title(r'$dB_Y/dt, \qquad \Sigma_{{total}} = {0:d},\; \Sigma_{{nz}} = {1:d}$'.format(sum(y_count),np.sum(np.array(y_count)>0)))
        gl = axes[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
    elif i == 2:
        cax = axes[i].scatter(z_phi,90.-np.array(z_theta), s = marker_size, c=z_count,cmap=cmap, vmin=0-0.5, vmax=(max_count)+0.5)
        axes[i].set_title(r'$dB_Z/dt, \qquad \Sigma_{{total}} = {0:d},\; \Sigma_{{nz}} = {1:d}$'.format(sum(z_count),np.sum(np.array(z_count)>0)))
        gl = axes[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        
    axes[i].coastlines()
    #gl = axes[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=[1,0,0,1],
    #                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-180, -135., -90, -45., 0, 45., 90, 135., 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

cb_ax = f.add_axes([0.15, -0.0, 0.7, 0.4])
#cb_ax.scatter([0,1,2,3],[0,1,2,3],[0,1,2,4])
cbar = f.colorbar(cax, ax=cb_ax, orientation = 'horizontal')
cb_ax.set_axis_off()

cbar.set_ticks(range(0,max_count+1))
outfname = 'Jerk'+str(jerk_number+1)+'_' + str(int(10*height_threshold)) + '_' + str(distance_threshold)
f.savefig(outfname + '.pdf',bbox_inches = 'tight')
f.savefig(outfname + '.png',bbox_inches = 'tight')
plt.close()

#Plot against timeseries at particular location

import chaosmagpy as cp
# import the dataset
import h5py
filepath = 'Gauss_Bsurf.mat'
arrays = {}
f = h5py.File(filepath,'r')
for k, v in f.items():
    arrays[k] = np.array(v)

coeffs = arrays['gnm'][:,:].T
time = arrays['timers'].flatten()
print( 'Shape of gmn array: ', arrays['gnm'].shape )
print ('Times (in years) of output', time )
radius = 6371.2
TIMES = time_yearly
NUM_DATA = len(TIMES)
TIMES_MIN = TIMES.min()
TIMES_MAX = TIMES.max()
CP_NBINS = 1*np.int(TIMES_MAX - TIMES_MIN) #one per year

plt.figure()
f, (ax1, ax2,ax3,ax4) = plt.subplots(4, 1, figsize=(13,9) )

for j in range(len(results)):
    theta = results[j][0]
    phi = results[j][1]
    component = results[j][2]
    CP = results[j][3]
    if(theta==90 and phi == 80):
        theta=90
        phi=80
        Br, Btheta, Bphi = cp.model_utils.synth_values(coeffs, radius, theta, phi,nmax=13)
        Br_yearly,Btheta_yearly, Bphi_yearly = np.interp(time_yearly, time, Br ), np.interp(time_yearly, time, Btheta ), np.interp(time_yearly, time, Bphi )
        Bx_dot, By_dot, Bz_dot = -np.gradient(Btheta_yearly,time_yearly), np.gradient(Bphi_yearly,time_yearly), -np.gradient(Br_yearly,time_yearly)
        left_edges = np.linspace(TIMES_MIN, TIMES_MAX, CP_NBINS,endpoint=False)
        if component == 0:
            ax1.bar(left_edges, CP, align='edge', width = 0.85*(left_edges[1] - left_edges[0]))
            ax1.set_xlim(time_yearly.min(), time_yearly.max() )
            ax5 = ax1.twinx()
            ax5.plot(TIMES,Bx_dot,'b')
            ax1.set_title(r'$d{B_X}/dt$')

        if component == 1:
            ax2.bar(left_edges, CP, align='edge', width = 0.85*(left_edges[1] - left_edges[0]))
            ax2.set_xlim(time_yearly.min(), time_yearly.max() )
            ax5 = ax2.twinx()
            ax5.plot(TIMES,By_dot,'b')
            ax2.set_title(r'$d{B_Y}/dt$')

        if component == 2:
            ax3.bar(left_edges, CP, align='edge', width = 0.85*(left_edges[1] - left_edges[0]))
            ax3.set_xlim(time_yearly.min(), time_yearly.max() )
            ax5 = ax3.twinx()
            ax5.plot(TIMES,Bz_dot,'b')
            ax3.set_title(r'$d{B_Z}/dt$')

        time_EJ, EJ = np.loadtxt('Jerk_energy.dat',unpack=True)
        ax4.plot(time_EJ, EJ )
        ax4.set_xlim(8860, 8900 )
        ax4.set_title('Jerk energy')


        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)

outfname = 'Jerk'+str(jerk_number+1)+'_ts_theta90_phi80'
plt.savefig(outfname+'.pdf')
plt.savefig(outfname+'.png')
