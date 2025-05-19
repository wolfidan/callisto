#%%
'''
This script process the solar scans automatically. It also performs a 2D paraboloid
fit in order to get the location of the Sun in the image and the hpbw.

Author: Andrea Francesco Battaglia

History:
  - 2025/03/27 [Andrea F. Battaglia]: Created.
  - 2025/04/02 [Andrea F. Battaglia]: Script adapted for the 6 solar scans (one every hour).
  - 2025/05/14 [Daniel Wolfensberger]: Complete rewrite to support arbitrary number of solar scans at any time
'''

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
#import time
import configparser
import traceback  # for tracing errors
import argparse
import datetime
import glob

import warnings
warnings.filterwarnings('ignore')

from astropy.io import fits
from scipy.optimize import curve_fit

from utils import read_motor_log, paraboloid_2d, gaussian_2d
from utils import convert_decimal_hour_timedelta, readfit, log_execution
from utils import get_fit_files, average_around_time

### Path where to store the diagnostics images
folder_images = 'C:\\xrt\\output\\solar_scan\\maps_solar-scan'

### Directlry where to store the pointing offsets
directory_csv_offsets = 'C:\\xrt\\output\\solar_scan'

### Path to the configsun file, where we get the aziref and eleref values
path_configsun = 'C:\\xrt\\src\\PythonScripts\\TrackingSun\\configsun.ini'

### Frequency to consider for the analysis
freq = 10637 # MHz

### Name of entry in motor log files that specifies sun scans
message_sunscans = " New scanning59 position"

### If two sunscans are further away in time than this time (in seconds) 
# we consider them as separate sunscans
max_time_diff_sunscans = 60

### If the distance of the center of the fitted paraboloid is larger than this value,
# we consider the measure as "off-pointing"
dist_quality_check = 0.7 # deg

# Log filename
logfname = 'C:\\xrt\\output\\solar_scan\\log_process_solar_scan.txt'

###########################################################################

day = datetime.datetime.now().strftime('%Y%m%d')
day_dt = datetime.datetime.strptime(day, "%Y%m%d")
fit_files, times_fit_files = get_fit_files(day)

### Read configuration file and get aziref and eleref
config = configparser.ConfigParser()
config.read(path_configsun)
aziref    = config.getfloat('Tracker','aziref')
eleref    = config.getfloat('Tracker','eleref')

### read txt file
df_motor = read_motor_log(day, path_log_motor)

# Get all scans from motor log
all_scans = df_motor[df_motor['Message'] == message_sunscans]
# Find where there are time breaks in scan entries
time_breaks = np.where(np.diff(all_scans['TIME [UT]']) > max_time_diff_sunscans/3600)[0]
# add zero
time_breaks = [-1] + list(time_breaks) + [len(all_scans) - 1]
nscans = len(time_breaks) - 1

log_execution(logfname, f"Found {nscans} scans in the motor log file for day {day}")

plt.rcParams.update({'font.size': 18})
nc = 2 # nb of plot columns
nr = int(np.ceil(nscans / nc))
fig1, ax1 = plt.subplots(nr,nc, figsize=(25, 7 * nr))
nc = 3 # nb of plot columns
nr = int(np.ceil(nscans / nc))
fig2, ax2 = plt.subplots(nr,nc, figsize=(26, 8 * nr), sharex=True, sharey=True)
ax1 = ax1.ravel()
ax2 = ax2.ravel()

# To be filled for every scan
all_quality_checks = []
all_params = []
all_start_times = []
all_end_times = []

for i in range(nscans): # loop on all sun raster scans
    # Get times of scan
    scan_data = all_scans[(time_breaks[i] + 1): (time_breaks[i+1]+1)]
    times = np.array(scan_data["TIME [UT]"])
    timedelta = convert_decimal_hour_timedelta(times)
    dtimes = [day_dt + dt for dt in timedelta]
    # add 1.1 second to every time to get time at middle of scan
    dtimes_with_offset = [dt + datetime.timedelta(seconds = 1.1) for dt in dtimes]
    log_execution(logfname, f"Processing scan from {dtimes[0]} to {dtimes[-1]}")
    try:
        # Find fit file closest in time (just before)
        idx_closest = np.where(np.array(times_fit_files) < dtimes[-1])[0][-1]
        
        # Read fit file
        fit_data = readfit(fit_files[idx_closest])
        log_execution(logfname, f"Corresponding fit file: {fit_files[idx_closest]}")
        
        ### Get the index of the frequency to look at
        idx_freq_closest = np.argmin(np.abs(fit_data["freqs"] - freq))

        ##### only select times relevant for the scan in the FIT files
        buffer = datetime.timedelta(seconds = 1.5)
        mask_inscan = np.logical_and(fit_data["datetimes"] >= dtimes[0] - buffer,
                            fit_data["datetimes"] <= dtimes[-1] + buffer)
        
        data_avg = average_around_time(fit_data["datetimes"],
                                    fit_data["data"][idx_freq_closest],
                                    dtimes_with_offset)
        idx_max = np.argmax(data_avg)
        
        # get relative azi and ele
        delta_azi = (scan_data["Azimuth[deg]"] - 
                    np.mean(scan_data["Azimuth[deg]"])).to_numpy()
        delta_ele = (scan_data["Elevation[deg]"] - 
                    np.mean(scan_data["Elevation[deg]"])).to_numpy()

        # correct the delta azimuth for the solid angle effects due to the elevation
        delta_azi = delta_azi * np.cos(np.deg2rad(scan_data["Elevation[deg]"].to_numpy()))
            
        # This is to define the ranges of the plot, plus some other stuff
        vmax = np.max([np.abs(delta_azi), np.abs(delta_ele)])
        vmax += 0.2
        vmin = -vmax
        
        # Do the fit
        ### Create meshgrid for smooth plotting
        x = np.linspace(-vmax, vmax, 100)
        y = np.linspace(-vmax, vmax, 100)
        X, Y = np.meshgrid(x, y)

        plot_fit = True
        quality_check = 1
        try:
            params, _ = curve_fit(paraboloid_2d, (delta_azi, delta_ele), data_avg, 
                p0=[np.max(data_avg), 0, 0, 1.5, 1.5],
                bounds=([0.01, -3, -3, 0.1, 0.1],
                [2*np.max(data_avg), 3, 3, 6, 6]))
            
            Z = paraboloid_2d((X, Y), *params)
            Z_max = np.max(Z)
            labels = {
                Z_max*0.9999: 'max',
                Z_max-3: 'hpbw'
            }
            dist = np.sqrt(params[1]**2+params[2]**2)
            quality_check = dist < dist_quality_check
        except Exception as err:
            plot_fit = False
            params = [np.nan, np.nan, np.nan, np.nan, np.nan]
            quality_check = 0
            log_execution(logfname, f"The fit of the scan at {dtimes[-1].strftime("%Y-%m-%d %H:%M")} has failed!")
            log_execution(logfname, f"Error is {err}")
            
        
        fit_data_in_scan = fit_data["data"][idx_freq_closest, mask_inscan]
        fit_times_in_scan = fit_data["datetimes"][mask_inscan]
        
        # First diagnostics plot: time profiles with extracted values

        dt_scan = 0.75
        fig1.suptitle(f"Solar scan at {freq} MHz on {day_dt.strftime("%Y-%m-%d")} (aziref = {aziref} deg; eleref = {eleref} deg)", 
                    y=0.98)

        ax1[i].plot(fit_times_in_scan, fit_data_in_scan, label='Morning', lw=2)
        ax1[i].scatter(dtimes_with_offset, data_avg, c='r', s=100, label='Morning (averaged)')
        ax1[i].set_ylabel('dBadu')
        ax1[i].set_xlabel('UT')
        ax1[i].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
        ax1[i].grid()
        ax1[i].tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax1[i].tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax1[i].legend()
        for j in range(len(dtimes_with_offset)):
            ax1[i].vlines(dtimes_with_offset[j], 0, 255, color='k', linestyle='dotted', alpha=0.5)
            ax1[i].axvspan(dtimes_with_offset[j]-pd.to_timedelta(dt_scan/2, unit='s'), 
                    dtimes_with_offset[j]+pd.to_timedelta(dt_scan/2, unit='s'), color='gray', alpha=0.4)
        ax1[i].set_ylim([np.min(fit_data_in_scan)-3, np.max(fit_data_in_scan)+3])
        ax1[i].xaxis.set_minor_locator(plt.matplotlib.dates.SecondLocator())
        ax1[i].tick_params(axis='x', which='minor', direction='in', 
                        top=True, bottom=True, width=1, length=4)

        # Second diagnostics plot: map of the scans
        fig2.suptitle(f"Solar scan: {dtimes[-1].strftime("%Y/%m/%d")} (aziref = {aziref} deg; eleref = {eleref} deg)",
                    y=0.98, fontsize=22, bbox=dict(facecolor='gray', alpha=0.2))

        ax2[i].plot(delta_azi, delta_ele, 'k--', alpha=0.4, lw=2)
        if plot_fit:
            contours_m_before = ax2[i].contour(X, Y, Z, 
                                    levels=[Z_max-3, Z_max*0.9999], 
                                    colors=['k', 'k'], alpha=0.5, linewidths=2)
            def custom_fmt(x):
                return labels.get(x, '')
            ax2[i].clabel(contours_m_before, inline=True, fontsize=16, fmt=custom_fmt)
            ax2[i].text(params[1]+0.1, params[2]+0.1, '({:.2f}, {:.2f})'.format(params[1], params[2]), fontsize=18, color='gray')
        cf_m_before = ax2[i].scatter(delta_azi, delta_ele, c=data_avg, s=100, cmap='hot', 
                        edgecolor='k', linewidth=1)
        ax2[i].set_xlabel(r'$\delta_{azimuth}$ [deg]')
        ax2[i].set_ylabel(r'$\delta_{elevation}$ [deg]')
        ax2[i].set_title(f"{dtimes[0].strftime("%H:%M:%S")}-{dtimes[-1].strftime("%H:%M:%S")} UT")
        ax2[i].set_xlim([-vmax, vmax])
        ax2[i].set_ylim([-vmax, vmax])
        ax2[i].set_xticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax2[i].set_yticks(np.arange(-vmax, vmax+0.1, 0.4))
        ax2[i].tick_params(axis='x', direction='in', top=True, bottom=True, width=1., length=7)
        ax2[i].tick_params(axis='y', direction='in', left=True, right=True, width=1., length=7)
        ax2[i].grid()
        cbar_m_before = plt.colorbar(cf_m_before)
        cbar_m_before.set_label('dBadu', labelpad=-59)
        ax2[i].add_patch(plt.Circle((delta_azi[idx_max], delta_ele[idx_max]), 0.12, color='r', fill=False, lw=2))
        ax2[i].text(delta_azi[idx_max], delta_ele[idx_max]+0.15, 'Max', color='r', fontsize=14, ha='center')
        if not quality_check:
            ax2[i].text(0, 0.3, 'off-pointing', fontsize=40, color='red', fontweight='bold', alpha=0.5, ha='center', va='center')
            
        fig1.subplots_adjust(top=0.93)
        fig2.subplots_adjust(top=0.93)
        # delete empty axes if needed
        if nscans % 2:
            ax1[-1].set_axis_off()
            ax2[-1].set_axis_off()
        
        #fig1.savefig(os.path.join(folder_images,f'solar-scan59_time-profiles_{day_dt.strftime("%Y-%m-%d")}.png'))
         #fig2.savefig(os.path.join(folder_images,f'solar-scan59_map_{day_dt.strftime("%Y-%m-%d")}.png'))

        # If everything runs successfully
        log_execution(logfname, "Successfully run")

    except Exception as err:
        ### Save NaNs in the csv file if it fails
        params = [np.nan, np.nan, np.nan, np.nan, np.nan]
        quality_check = 0
        
        ### If an error occurs, log the error message and stack trace
        error_message = f"Error: {str(err)}\n{traceback.format_exc()}"
        log_execution(logfname, error_message)
    finally:
        all_params.append(params)
        all_quality_checks.append(quality_check)
        all_start_times.append(dtimes[0].strftime("%Y-%m-%d %H:%M:%S"))
        all_end_times.append(dtimes[-1].strftime("%Y-%m-%d %H:%M:%S"))
        
### Find the csv with the pointing offsets
YYYY = all_start_times[0][0:4]
filename_csv_offsets = 'pointing-offsets_%s.csv' % YYYY
filepath_csv_offsets = os.path.join(directory_csv_offsets, filename_csv_offsets)
# if it does not exist, create it
if os.path.exists(filepath_csv_offsets):
    df_offsets = pd.read_csv(filepath_csv_offsets)
else:
    df_offsets = pd.DataFrame(columns=['start_scan_UTC', 'end_scan_UTC', 'amplitude_fit', 'delta_azi_fit', 'delta_ele_fit', 'sigma_azi_fit', 'sigma_ele_fit', 'aziref', 'eleref', 'quality_check'])
    df_offsets.to_csv(filepath_csv_offsets, index=False)

### Store the results in the update the csv file in df_flux by adding the new values
df_offsets = pd.read_csv(filepath_csv_offsets)
for i in range(len(all_start_times)):
    new_row = pd.DataFrame({'start_scan_UTC': str(all_start_times[i]), 
        'end_scan_UTC': str(all_end_times[i]), 
        'amplitude_fit': all_params[i][0], 
        'delta_azi_fit': all_params[i][1],
        'delta_ele_fit': all_params[i][2], 
        'sigma_azi_fit': all_params[i][3], 
        'sigma_ele_fit': all_params[i][4],
        'aziref': aziref,
        'eleref': eleref,
        'quality_check': all_quality_checks[i]}, index=[0])

    df_offsets = pd.concat([df_offsets,new_row], ignore_index=True)
df_offsets = df_offsets.sort_values(by='start_scan_UTC')
df_offsets = df_offsets.drop_duplicates(subset="start_scan_UTC")
df_offsets.to_csv(filepath_csv_offsets, index=False)


