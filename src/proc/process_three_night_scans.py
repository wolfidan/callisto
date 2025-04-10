# -*- coding: utf-8 -*-

'''
This script automatically process the three night scans that are run in sequence
  - from 21:15 UT to 23:15 UT -> TV-sat scan
  - from 23:30 UT to 23:55 UT -> zenith-ground scan
  - from 00:00 UT to 04:23 UT -> full-sky scan
Summary images will be produced as output.

It also includes useful functions to properly read the scan data.
These functions distinguish between the different types of scan.
There are also some functions to allow for the plotting of context (e.g., the 
geostationary orbit, the mountains, ...) in the plot of the scan data.

Author: Andrea F. Battaglia

History:
  - 2025/03/25: Created.
'''


import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import math
import ephem
import traceback  # for tracing errors
import datetime
import argparse
import glob

from astropy.io import fits
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

import warnings
warnings.filterwarnings('ignore')

### Location parameter for Locarno Monti
MyHome = ephem.Observer()
MyHome.lon      = '8.7875'    # east +°
MyHome.lat      = '46.1723'  # north +°
MyHome.elev     = 374       # altitude in m asl
MyHome.temp     = 15        # °C
MyHome.pressure = 900       # mbar
name            = 'Locarno Monti' # name of the location from where the coordinates are taken


#**********************************************************************

def log_execution(message, plot_separation=True):
    """
    Description:
        This is to store in an external file a log, whether the script run
        successfully or an error occurred.
    """
    with open('C:\\xrt\\output\\night_scan\\log_process_night_scan.txt', 'a') as log_file:
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_file.write(f'\n{timestamp} - {message}')
        if plot_separation == True: 
            log_file.write(f'\n************************************************\n')
        
#**********************************************************************

def readfit(path):
    
    with fits.open(path) as hdu:
        # https://docs.astropy.org/en/stable/io/fits/
        dict_fitfile = {
            'data'       : hdu[0].data, # .astype(np.uint8)
            'freqs'      : hdu[1].data  ['Frequency'][0], # extract frequency axis
            'timeax'     : hdu[1].data  ['Time'][0],      # extract time axis
            'dT'         : hdu[0].header['CDELT1'],       # extract time resolution
            'date'       : hdu[0].header['DATE-OBS'],     # take first file
            'T0'         : hdu[0].header['TIME-OBS'],     # take first file
            'instrument' : hdu[0].header['INSTRUME'],
            'content'    : hdu[0].header['CONTENT']}
        
    return dict_fitfile

#**********************************************************************

def T0fits(files):
    nfits = len(files)
    T0 = []
    for i in range(nfits):
        FITfile = readfit(files[i])
        T0.append(FITfile['T0'])
    return T0

#**********************************************************************

def readcsv_scan(path_csv, date_scan):
    
    df_scan = pd.read_csv(path_csv)
    time_scan = df_scan['time_UT']
    time_scan = [pd.to_datetime(date_scan+' '+t) for t in time_scan] 
    azi_scan = df_scan['azimuth_deg']
    ele_scan = df_scan['elevation_deg']
    
    return time_scan, azi_scan, ele_scan

#**********************************************************************

def readscan(files, time_scan):
    '''
    This function assumes that files contains the list of FIT files of the full sky scan.
    No files associated with other scans have to be present in the list.

    Parameters:
    ----------
        - files: list of FIT files associated with the full-sky scan
        - time_scan: list of time stamps of the scan (from the csv file)

    Returns:
    --------
        - data_scan: 2D array [freq, time] with the data of the scan (in ADC units)
        - freqs_MHz: 1D array with the frequency axis
        - all_times: all times of the FIT file 
        - all_data:  all data of the FIT file
    '''
    
    nfiles = len(files)

    ### open the first file to get the time and freq axes
    file0 = readfit(files[0])
    freqs_MHz = file0['freqs']
    timeax = file0['timeax']
    nfreq = len(freqs_MHz)
    ntimes = len(timeax)

    ### Get data and times of all time stamps in all files
    all_data = np.zeros((nfreq, nfiles*ntimes-20))
    all_times = []
    dates = [pd.to_datetime(readfit(file)['date'].replace('/', '-')) for file in files]
    for i, (file, base_date) in enumerate(zip(files, dates)):
        FITfile = readfit(file)
        timeax_FIT = FITfile['timeax']
        T0_FIT = pd.to_datetime(FITfile['T0']).time()
        
        base_datetime = pd.Timestamp.combine(base_date.date(), T0_FIT)
        time_axis_FIT = base_datetime + pd.to_timedelta(timeax_FIT, unit='s')
        
        if i != len(files) - 1:
            all_data[:, i*ntimes:(i+1)*ntimes] = FITfile['data']
            all_times.extend(time_axis_FIT)
        else:
            all_data[:, i*ntimes:(i+1)*ntimes-20] = FITfile['data'][:,0:3600-20]
            all_times.extend(time_axis_FIT[0:3600-20])
    
    ##### In the old version, we converted everything in dB, now we leave it as ADC units    
    ### conversion to dB
    #all_dB = all_data/255.*2500.0/25.4

    ### Average data over time stamps between imin and imax
    data_scan = np.zeros((nfreq, len(time_scan)))
    imin = 40
    imax = 200
    time_scan_np = np.array(time_scan).astype('datetime64[ns]').astype('float64')
    all_times_np = np.array(all_times).astype('datetime64[ns]').astype('float64')
    differences = np.abs(time_scan_np.reshape(-1, 1) - all_times_np)
    indices = np.argmin(differences, axis=1)
    for i, idx in enumerate(indices):
        data_scan[:,i] = all_data[:,idx+imin:idx+imax].mean(axis=1)

    return data_scan, freqs_MHz, all_times, all_data

#**********************************************************************

def GeoSat(sat_lon):
    lat, lon = math.degrees(float(MyHome.lat)), math.degrees(float(MyHome.lon))
    rlat = math.radians(lat)
    rlon = math.radians(lon)
    rsat = math.radians(sat_lon)
    L = rsat - rlon
    D = math.acos(math.cos(rlat) * math.cos(L))
    az = math.degrees(math.acos(-math.tan(rlat) / math.tan(D)))
    az = az if L > 0 else 360 - az
    cd = math.cos(D)
    num = cd - 1 / 6.62
    den = (1 - cd * cd) ** 0.5
    el = math.degrees(math.atan(num / den))
    return(az,el)

#**********************************************************************

def plot_context(ax):
    '''
    This function sets the context for the plot of the scan data (in the case of
    2D plot with azimuth on the x-axis and elevation on the y-axis).
    
    Parameters:
    ----------
        - ax: axis of the plot
    '''
    
    path_csv_mountains = 'C:\\xrt\\src\\proc\\profiles_mountains_LocarnoMonti.csv'

    ### read csv file with the profile of the mountains
    df_mountains = pd.read_csv(path_csv_mountains, sep=',')
    azi_mountains = df_mountains['azimuth_deg'].values
    height_mountains = df_mountains['height_m'].values - 373 # 373 is the altitude of Locarno Monti
    distance_mountains = df_mountains['distance_km'].values * 1000
    elevation_mountains = np.arctan(height_mountains/distance_mountains) * 180 / np.pi

    ### Calculate the azimuth and elevation of sun at the summer and winter solstice
    year_ss = 2025
    month_ss = 6
    day_ss = 21

    AZI_ss = []
    ELE_ss = []
    HM_ss  = []
    hh_ss  = []
    hhmm_ss = []
    date_ss = []

    sun_ss = ephem.Sun()
    for h in range(0,24):
        for m in range(0,60,10):
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year_ss,month_ss,day_ss,h,m,0)
            sun_ss.compute(MyHome)
            Sazi = math.degrees(sun_ss.az)
            Sele = math.degrees(sun_ss.alt)
            
            #plt.plot(Sazi,Sele,'.r')
            hm = h+m/60
            HM_ss = np.append(HM_ss,hm)
            AZI_ss = np.append(AZI_ss,Sazi)
            ELE_ss = np.append(ELE_ss,Sele)
            st = '{:4.0f}'.format(hm)
            hh_ss = np.append(hh_ss,st)
            hhmm_ss = np.append(hhmm_ss,'{:02d}:{:02d}'.format(h,m))
            date_ss = np.append(date_ss,'{:4d}-{:02d}-{:02d}'.format(year_ss,month_ss,day_ss))

    year_ws = 2025
    month_ws = 12
    day_ws = 21

    AZI_ws = []
    ELE_ws = []
    HM_ws  = []
    hh_ws  = []
    hhmm_ws = []
    date_ws = []

    sun_ws = ephem.Sun()
    for h in range(0,24):
        for m in range(0,60,10):
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year_ws,month_ws,day_ws,h,m,0)
            sun_ws.compute(MyHome)
            Sazi = math.degrees(sun_ws.az)
            Sele = math.degrees(sun_ws.alt)
            
            #plt.plot(Sazi,Sele,'.r')
            hm = h+m/60
            HM_ws = np.append(HM_ws,hm)
            AZI_ws = np.append(AZI_ws,Sazi)
            ELE_ws = np.append(ELE_ws,Sele)
            st = '{:4.0f}'.format(hm)
            hh_ws = np.append(hh_ws,st)
            hhmm_ws = np.append(hhmm_ws,'{:02d}:{:02d}'.format(h,m))
            date_ws = np.append(date_ws,'{:4d}-{:02d}-{:02d}'.format(year_ws,month_ws,day_ws))

    az_sat = []
    el_sat = []
    for sat in range (0,360,5):
        a,e = GeoSat(sat)
        if (e > 0):
                ax.plot(a,e,'.',c='gray')
                az_sat = np.append(az_sat,a)
                el_sat = np.append(el_sat,e)
        ax.plot(360-a,e,'.',c='gray') 
        az_sat = np.append(az_sat,360-a)
        el_sat = np.append(el_sat,e)
    ax.plot(a,e,'.',c='gray',label='Geostationary orbit')
    ax.plot(AZI_ss,ELE_ss,'--',c='black',alpha=0.8)
    ax.plot(AZI_ws,ELE_ws,'--',c='black',alpha=0.8)
    ax.plot(azi_mountains,elevation_mountains,'-',c='black')

#**********************************************************************

def calibrate_night_scan(data_scan_TVsat, data_scan_zg, data_scan_fullsky, ele_scan_zg, freqs_MHz):
    '''
    This function should calibrate the night scan observations.
    HOWEVER, IT DOES NOT SEEM TO GIVE MEANINGFUL RESULTS!!!! This could be
    due to the assumption we do for the calibration.

    AT THE MOMENT IT IS NOT USED!
    '''

    ### Elevation angles for the calibration Thot and Tcold
    elevation_Tcold = 77
    elevation_Thot = -14
    ### for now, we just use a hard coded value for Thot and Tcold, but we should use the values as
    ## we use for the daily flux solar calibration
    Thot_K = 273
    Tcold_K = 20

    idx_Tcold = np.argmin(np.abs(ele_scan_zg - elevation_Tcold))
    idx_Thot = np.argmin(np.abs(ele_scan_zg - elevation_Thot))
    #Ihot = np.mean(data_scan_zg[:, idx_Thot], axis=0)
    #Icold = np.mean(data_scan_zg[:, idx_Tcold], axis=0)
    Ihot = data_scan_zg[:, idx_Thot]   ## this is in ADC units
    Icold = data_scan_zg[:, idx_Tcold]  ## this is in ADC units

    ### Convert all data into dB
    data_scan_TVsat_dB = data_scan_TVsat/255.*2500.0/25.4
    data_scan_zg_dB = data_scan_zg/255.*2500.0/25.4
    data_scan_fullsky_dB = data_scan_fullsky/255.*2500.0/25.4
    Ihot_dB = Ihot/255.*2500.0/25.4
    Icold_dB = Icold/255.*2500.0/25.4

    ### convert to linear
    data_scan_TVsat_lin = 10**(data_scan_TVsat_dB/10.)
    data_scan_zg_lin = 10**(data_scan_zg_dB/10.)
    data_scan_fullsky_lin = 10**(data_scan_fullsky_dB/10.)
    Ihot_lin = 10**(Ihot_dB/10.)
    Icold_lin = 10**(Icold_dB/10.)

    ### Simplified version: we assume that the ground and the mountains are far away, so we can use the usual formula
    T_TVsat = np.zeros((len(freqs_MHz), len(time_scan_TVsat)))
    for tt in range(len(data_scan_TVsat_lin[0,:])):
        T_TVsat[:, tt] = (Thot_K - Tcold_K) * (data_scan_TVsat_lin[:,tt] - Icold_lin) / (Ihot_lin - Icold_lin)

    T_zg = np.zeros((len(freqs_MHz), len(ele_scan_zg)))
    for tt in range(len(data_scan_zg_lin[0,:])):
        T_zg[:, tt] = (Thot_K - Tcold_K) * (data_scan_zg_lin[:, tt] - Icold_lin) / (Ihot_lin - Icold_lin)

    T_fullsky = np.zeros((len(freqs_MHz), len(ele_scan_fullsky)))
    for tt in range(len(data_scan_fullsky_lin[0,:])):
        T_fullsky[:, tt] = (Thot_K - Tcold_K) * (data_scan_fullsky_lin[:,tt] - Icold_lin) / (Ihot_lin - Icold_lin)

    return T_TVsat, T_zg, T_fullsky

#**********************************************************************

def gaussian_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    x, y = coords
    return amplitude * np.exp(-(((x-x0)**2)/(2*sigma_x**2) + ((y-y0)**2)/(2*sigma_y**2)))

#**********************************************************************



def main():

    print()
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('Processing of the night observations. It takes about 10 seconds.')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print()

    ### Path to the csv with the times and coordinates of the scans
    path_csv_TVscan = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_TV-sat.csv'
    path_csv_zg = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_zenith-ground.csv'
    path_csv_fullsky = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_full-sky.csv'

    ### Path to the folder containing the FIT files of the night-sky scan
    path_folder_FIT = 'C:\\xrt\\output\\data\\raw\\FITfiles'

    ### Where to store the plots
    folder_store = 'C:\\xrt\\output\\night_scan\\scan-results\\'

    ### Frequencies to check the night-sky scan [MHz]
    in_freq = [10640, 10800, 10965, 11065]
    
    ### We use this frequency to do the fit of the TV-sat
    freq_sat = 10797 # MHz

    try:
    
        ### The user can define a day in the past for re-analysis
        # By default, today is taken
        parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        epilog='')
        
        parser.add_argument(
            '-day', type=str, metavar=str.__name__, 
            action='store', required=False, default=None,
            help='Choose day to analyse (default: today), format: YYYYMMDD')
        
        args = parser.parse_args()
        day = args.day

        ### If a day was not specified, use today's date
        if not day: day = datetime.datetime.now().strftime('%Y%m%d')

        day_pandas = pd.to_datetime(day, format='%Y%m%d')
        day_pandas_yesterday = day_pandas - pd.Timedelta(days=1)
        yesterday = day_pandas_yesterday.strftime('%Y%m%d')
                
        filename_format_today = 'meteoswiss_%s_*_01.fit' % day #MyFile = 'SWISS-METEO_20231127_*_01.fit'
        filename_format_yesterday = 'meteoswiss_%s_*_01.fit' % yesterday #MyFile = 'SWISS-METEO_20231127_*_01.fit'
        filepath_format_today = os.path.join(path_folder_FIT, filename_format_today)
        filepath_format_yesterday = os.path.join(path_folder_FIT, filename_format_yesterday)
        filepaths_today = sorted(glob.glob(filepath_format_today))
        filepaths_yesterday = sorted(glob.glob(filepath_format_yesterday))
                    
        ### Get the starting time of each FIT file
        T0_fits_yesterday = T0fits(filepaths_yesterday)
        T0_fits_today = T0fits(filepaths_today)
        nfiles_yesterday = len(T0_fits_yesterday)
        nfiles_today = len(T0_fits_today)

        ### Sort out the files by time of the day
        TVsat_files = []
        zg_files = []
        fullsky_files = []
        for i in range(nfiles_yesterday):

            if T0_fits_yesterday[i] > '21:15' and T0_fits_yesterday[i] < '23:16':
                TVsat_files.append(filepaths_yesterday[i])

            if T0_fits_yesterday[i] > '23:29' and T0_fits_yesterday[i] < '23:59':
                zg_files.append(filepaths_yesterday[i])

        for i in range(nfiles_today):

            if T0_fits_today[i] >= '00:00' and T0_fits_today[i] < '04:30':
                fullsky_files.append(filepaths_today[i])

        TVsat_files.sort()
        zg_files.sort()
        fullsky_files.sort()

        ### Dates scan
        date_scan_TVsat = readfit(TVsat_files[0])['date'].replace('/', '-')
        date_scan_zg = readfit(zg_files[0])['date'].replace('/', '-')
        date_scan_fullsky = readfit(fullsky_files[0])['date'].replace('/', '-')

        ### Read the CSVs of the scans
        time_scan_TVsat, azi_scan_TVsat, ele_scan_TVsat = readcsv_scan(path_csv_TVscan, date_scan_TVsat)
        time_scan_zg, azi_scan_zg, ele_scan_zg = readcsv_scan(path_csv_zg, date_scan_zg)
        time_scan_fullsky, azi_scan_fullsky, ele_scan_fullsky = readcsv_scan(path_csv_fullsky, date_scan_fullsky)

        ### Read the FIT files at the times of the scan and average the values for 40 seconds 
        data_scan_TVsat, freqs_MHz, times_FIT_TVsat, data_FIT_TVsat = readscan(TVsat_files, time_scan_TVsat)
        data_scan_zg, freqs_MHz, times_FIT_zg, data_FIT_zg = readscan(zg_files, time_scan_zg)
        data_scan_fullsky, freqs_MHz, times_FIT_fullsky, data_FIT_fullsky = readscan(fullsky_files, time_scan_fullsky)

        ### Selected frequencies
        idx_freqs = [np.argmin(np.abs(freqs_MHz - f)) for f in in_freq]

        ### calibrate the observations to get the brightness temperature
        #T_TVsat, T_zg, T_fullsky = calibrate_night_scan(data_scan_TVsat, data_scan_zg, data_scan_fullsky, ele_scan_zg, freqs_MHz)

        ### Ranges of the colorbars 
        # (not actually used, because we want the same colorbar every day, to better see the changes)
        #ADU_max_TVsat = np.max([data_scan_TVsat[idx_freqs[0],:], data_scan_TVsat[idx_freqs[1],:], data_scan_TVsat[idx_freqs[2],:], data_scan_TVsat[idx_freqs[3],:]]) 
        #ADU_min_TVsat = np.min([data_scan_TVsat[idx_freqs[0],:], data_scan_TVsat[idx_freqs[1],:], data_scan_TVsat[idx_freqs[2],:], data_scan_TVsat[idx_freqs[3],:]])
        #ADU_max_zg = np.max([data_scan_zg[idx_freqs[0],:], data_scan_zg[idx_freqs[1],:], data_scan_zg[idx_freqs[2],:], data_scan_zg[idx_freqs[3],:]])
        #ADU_min_zg = np.min([data_scan_zg[idx_freqs[0],:], data_scan_zg[idx_freqs[1],:], data_scan_zg[idx_freqs[2],:], data_scan_zg[idx_freqs[3],:]])
        #ADU_max_fullsky = np.max([data_scan_fullsky[idx_freqs[0],:], data_scan_fullsky[idx_freqs[1],:], data_scan_fullsky[idx_freqs[2],:], data_scan_fullsky[idx_freqs[3],:]])
        #ADU_min_fullsky = np.min([data_scan_fullsky[idx_freqs[0],:], data_scan_fullsky[idx_freqs[1],:], data_scan_fullsky[idx_freqs[2],:], data_scan_fullsky[idx_freqs[3],:]])

        #K_max_TVsat = np.max([T_TVsat[idx_freqs[0],:], T_TVsat[idx_freqs[1],:], T_TVsat[idx_freqs[2],:], T_TVsat[idx_freqs[3],:]])
        #K_min_TVsat = np.min([T_TVsat[idx_freqs[0],:], T_TVsat[idx_freqs[1],:], T_TVsat[idx_freqs[2],:], T_TVsat[idx_freqs[3],:]])
        #K_max_zg = np.max([T_zg[idx_freqs[0],:], T_zg[idx_freqs[1],:], T_zg[idx_freqs[2],:], T_zg[idx_freqs[3],:]])
        #K_min_zg = np.min([T_zg[idx_freqs[0],:], T_zg[idx_freqs[1],:], T_zg[idx_freqs[2],:], T_zg[idx_freqs[3],:]])
        #K_max_fullsky = np.max([T_fullsky[idx_freqs[0],:], T_fullsky[idx_freqs[1],:], T_fullsky[idx_freqs[2],:], T_fullsky[idx_freqs[3],:]])
        #K_min_fullsky = np.min([T_fullsky[idx_freqs[0],:], T_fullsky[idx_freqs[1],:], T_fullsky[idx_freqs[2],:], T_fullsky[idx_freqs[3],:]])


        #**********************************************************************
        # Plot of the TV-sat (zoom-in version + gaussian fit)
        #**********************************************************************


        vmin = 138 #ADU_min_TVsat
        vmax = 185 #ADU_max_TVsat
        xlim = [155, 175]
        ylim = [26, 46]

        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(16, 14))

        fig.suptitle('CALLISTO TV-sat scanning mode (zoom-in) on '+date_scan_TVsat, fontsize=24, y=0.95)

        cbar_ax = fig.add_axes([0.11, 0.06, 0.78, 0.013])

        gs = fig.add_gridspec(2, 2, 
                            top=0.92,     # Higher top margin for title
                            bottom=0.12,  # Higher bottom margin for colorbar
                            left=0.1, 
                            right=0.9,
                            wspace=0.05, 
                            hspace=0.05)

        ax0 = fig.add_subplot(gs[0, 0])
        cf0 = ax0.scatter(azi_scan_TVsat, ele_scan_TVsat, c=data_scan_TVsat[idx_freqs[0],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=12)
        ax0.set_ylabel('Elevation [deg]')
        plot_context(ax0)
        ax0.set_xlim(xlim)
        ax0.set_ylim(ylim)
        ax0.text(156, 44.5, str(int(freqs_MHz[idx_freqs[0]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax0.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelbottom=False)
        #ax0.text(222,17,'Geostat. orbit',color='gray',fontsize=14,rotation=-45)
        #ax0.text(218,45,'Summer solstice',color='black',fontsize=14,rotation=-37)
        #ax0.text(200,7,'Winter solstice',color='black',fontsize=14,rotation=-37)
        #ax0.text(175,4,'Mountains',color='black',fontsize=14,rotation=-17)
        #ax0.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax1 = fig.add_subplot(gs[0, 1])
        cf1 = ax1.scatter(azi_scan_TVsat, ele_scan_TVsat, c=data_scan_TVsat[idx_freqs[1],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=12)
        plot_context(ax1)
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.text(156, 44.5, str(int(freqs_MHz[idx_freqs[1]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax1.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False, labelbottom=False)
        #ax1.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax2 = fig.add_subplot(gs[1, 0])
        cf2 = ax2.scatter(azi_scan_TVsat, ele_scan_TVsat, c=data_scan_TVsat[idx_freqs[2],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=12)
        ax2.set_xlabel('Azimuth [deg]')
        ax2.set_ylabel('Elevation [deg]')
        plot_context(ax2)
        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim)
        ax2.text(156, 44.5, str(int(freqs_MHz[idx_freqs[2]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax2.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True)
        #ax2.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        #ax2.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        ax3 = fig.add_subplot(gs[1, 1])
        cf3 = ax3.scatter(azi_scan_TVsat, ele_scan_TVsat, c=data_scan_TVsat[idx_freqs[3],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=12)
        ax3.set_xlabel('Azimuth [deg]')
        plot_context(ax3)
        ax3.set_xlim(xlim)
        ax3.set_ylim(ylim)
        ax3.text(156, 44.5, str(int(freqs_MHz[idx_freqs[3]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax3.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False)
        #ax3.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        #ax3.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        cbar = plt.colorbar(cf0, cax=cbar_ax, orientation='horizontal', ticklocation='bottom')
        cbar.set_label('ADU', y=-1.8)  # Move label up by adjusting y

        plt.tight_layout()
        plt.savefig(folder_store+'TV-sat-scan_zoom-in_'+str(date_scan_TVsat)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()


        ### Gauussian fit
        
        idx_TVsat = np.argmin(np.abs(freqs_MHz - freq_sat))
        data_scan_TVsat2fit = data_scan_TVsat[idx_TVsat,:]

        xlim = [160, 170]
        ylim = [31, 41]

        # Create meshgrid for smooth plotting
        x = np.linspace(xlim[0], xlim[1], 100)
        y = np.linspace(ylim[0], ylim[1], 100)
        X, Y = np.meshgrid(x, y)

        ### Do Gaussian fit
        params_fit, _ = curve_fit(gaussian_2d, (azi_scan_TVsat, ele_scan_TVsat), data_scan_TVsat2fit.flatten(), 
                                    p0=[np.max(data_scan_TVsat2fit), 165, 36, 0.5, 0.5])
        Z_fit = gaussian_2d((X, Y), *params_fit)
        Z_fit_norm = (Z_fit - np.min(Z_fit)) / (np.max(Z_fit) - np.min(Z_fit))
        max_fit = np.nanmax(Z_fit_norm)

        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(10, 8))

        fig.suptitle('CALLISTO TV-sat scanning mode with 2D Gaussian fit on '+date_scan_TVsat, fontsize=18, y=0.95)

        ax0 = plt.subplot(111)
        contours_n = ax0.contour(X, Y, Z_fit_norm, 
                                levels=[max_fit*0.5, max_fit*0.999], 
                                colors=['k', 'k', 'k'], alpha=0.5, linewidths=2)
        cf0 = ax0.scatter(azi_scan_TVsat, ele_scan_TVsat, c=data_scan_TVsat2fit, cmap='gnuplot', vmin=np.min(Z_fit), vmax=vmax, lw=8)
        ax0.clabel(contours_n, inline=True, fontsize=10, fmt='%.2f')
        ax0.set_ylabel('Elevation [deg]')
        ax0.set_xlabel('Azimuth [deg]')
        plot_context(ax0)
        ax0.set_xlim([xlim[0]-2, xlim[1]+2])
        ax0.set_ylim([ylim[0]-2, ylim[1]+2])
        ax0.text(158.4, 42.0, str(int(freqs_MHz[idx_TVsat]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax0.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True)
        ax0.grid(alpha=0.3)

        cbar = plt.colorbar(cf0)
        cbar.set_label('ADU') 
        ax0.text(164, 42, 'Center Gaussian: ('+str(round(params_fit[1], 1))+', '+str(round(params_fit[2], 1))+')', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        #ax0.text(167, 45, 'Semi-major axis: '+str(round(params_fit[3], 2)), fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        #ax0.text(167, 44, 'Semi-minor axis: '+str(round(params_fit[4], 2)), fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))

        plt.tight_layout()
        plt.savefig(folder_store+'TV-sat-scan_gaussian-fit_'+str(date_scan_TVsat)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()
        

        #**********************************************************************
        # Plot of the zenith-ground scan 
        #**********************************************************************


        ### MAP
        vmin = 138 #ADU_min_TVsat
        vmax = 185 #ADU_max_TVsat

        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(16, 14))

        fig.suptitle('CALLISTO zenith-ground scanning mode on '+date_scan_zg, fontsize=24, y=0.95)

        cbar_ax = fig.add_axes([0.11, 0.06, 0.78, 0.013])

        gs = fig.add_gridspec(2, 2, 
                            top=0.92,     # Higher top margin for title
                            bottom=0.12,  # Higher bottom margin for colorbar
                            left=0.1, 
                            right=0.9,
                            wspace=0.05, 
                            hspace=0.05)

        ax0 = fig.add_subplot(gs[0, 0])
        cf0 = ax0.scatter(azi_scan_zg, ele_scan_zg, c=data_scan_zg[idx_freqs[0],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=10)
        ax0.set_ylabel('Elevation [deg]')
        plot_context(ax0)
        ax0.set_xlim(105, 255)
        ax0.set_ylim(-5, 75)
        ax0.text(110, 70, str(int(freqs_MHz[idx_freqs[0]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax0.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelbottom=False)
        ax0.text(222,17,'Geostat. orbit',color='gray',fontsize=14,rotation=-45)
        ax0.text(218,45,'Summer solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(200,7,'Winter solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(175,4,'Mountains',color='black',fontsize=14,rotation=-17)
        ax0.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax1 = fig.add_subplot(gs[0, 1])
        cf1 = ax1.scatter(azi_scan_zg, ele_scan_zg, c=data_scan_zg[idx_freqs[1],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=10)
        plot_context(ax1)
        ax1.set_xlim(105, 255)
        ax1.set_ylim(-5, 75)
        ax1.text(110, 70, str(int(freqs_MHz[idx_freqs[1]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax1.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False, labelbottom=False)
        ax1.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax2 = fig.add_subplot(gs[1, 0])
        cf2 = ax2.scatter(azi_scan_zg, ele_scan_zg, c=data_scan_zg[idx_freqs[2],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=10)
        ax2.set_xlabel('Azimuth [deg]')
        ax2.set_ylabel('Elevation [deg]')
        plot_context(ax2)
        ax2.set_xlim(105, 255)
        ax2.set_ylim(-5, 75)
        ax2.text(110, 70, str(int(freqs_MHz[idx_freqs[2]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax2.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True)
        ax2.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax2.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        ax3 = fig.add_subplot(gs[1, 1])
        cf3 = ax3.scatter(azi_scan_zg, ele_scan_zg, c=data_scan_zg[idx_freqs[3],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=10)
        ax3.set_xlabel('Azimuth [deg]')
        plot_context(ax3)
        ax3.set_xlim(105, 255)
        ax3.set_ylim(-5, 75)
        ax3.text(110, 70, str(int(freqs_MHz[idx_freqs[3]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax3.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False)
        ax3.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax3.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        cbar = plt.colorbar(cf0, cax=cbar_ax, orientation='horizontal', ticklocation='bottom')
        cbar.set_label('ADU', y=-1.8)  # Move label up by adjusting y

        plt.tight_layout()
        plt.savefig(folder_store+'zenith-ground-scan_'+str(date_scan_zg)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()


        ### ELEVATION
        plt.rcParams.update({'font.size': 16})
        fig = plt.figure(figsize=(10, 8))
        fig.suptitle('CALLISTO zenith-ground scanning mode on '+date_scan_zg, y=0.95)

        ax = fig.add_subplot(111)
        ax.scatter(ele_scan_zg, data_scan_zg[idx_freqs[0],:], label=str(int(freqs_MHz[idx_freqs[0]]))+' MHz')
        ax.scatter(ele_scan_zg, data_scan_zg[idx_freqs[1],:], label=str(int(freqs_MHz[idx_freqs[1]]))+' MHz')
        ax.scatter(ele_scan_zg, data_scan_zg[idx_freqs[2],:], label=str(int(freqs_MHz[idx_freqs[2]]))+' MHz')
        ax.scatter(ele_scan_zg, data_scan_zg[idx_freqs[3],:], label=str(int(freqs_MHz[idx_freqs[3]]))+' MHz')
        ax.tick_params(axis='both', which='major', labelsize=16, direction='in',
                        length=5, width=2, colors='black', top=True, right=True)
        ax.set_xticks([-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
        ax.set_xticklabels(['-20', '-10', '0', '10', '20', '30', '40', '50', '60', '70', '80', '90'])
        ax.grid(alpha=0.3)

        ax.fill_between([30.1, 30.9], 100-15, 200+15, color='gray', alpha=0.4, label='Geostat. orbit')
        ax.fill_between([-15, 9], 100-15, 200+15, color='green', alpha=0.1, label='Ground')
        ax.fill_between([76.5, 78.5], 100-15, 200+15, color='blue', alpha=0.2, label='Tcold')
        ax.fill_between([-15, -13], 100-15, 200+15, color='red', alpha=0.2, label='Thot')

        ax.set_xlabel('Elevation [deg]')
        ax.set_ylabel('ADU')
        ax.set_ylim(138, 160)
        ax.set_xlim(-20, 95)
        ax.legend()

        plt.tight_layout()
        plt.savefig(folder_store+'zenith-ground-scan_elevation_'+str(date_scan_zg)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()

        
        #**********************************************************************
        # Plot of the full-sky scan
        #**********************************************************************


        ### OBSERVED
        vmin = 138 #ADU_min_TVsat
        vmax = 185 #ADU_max_TVsat

        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(16, 14))

        fig.suptitle('CALLISTO full-sky scanning mode on '+date_scan_fullsky, fontsize=24, y=0.95)

        cbar_ax = fig.add_axes([0.11, 0.06, 0.78, 0.013])

        gs = fig.add_gridspec(2, 2, 
                            top=0.92,     # Higher top margin for title
                            bottom=0.12,  # Higher bottom margin for colorbar
                            left=0.1, 
                            right=0.9,
                            wspace=0.05, 
                            hspace=0.05)

        ax0 = fig.add_subplot(gs[0, 0])
        cf0 = ax0.scatter(azi_scan_fullsky, ele_scan_fullsky, c=data_scan_fullsky[idx_freqs[0],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=8)
        ax0.set_ylabel('Elevation [deg]')
        plot_context(ax0)
        ax0.set_xlim(105, 255)
        ax0.set_ylim(-5, 75)
        ax0.text(110, 70, str(int(freqs_MHz[idx_freqs[0]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax0.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelbottom=False)
        ax0.text(222,17,'Geostat. orbit',color='gray',fontsize=14,rotation=-45)
        ax0.text(218,45,'Summer solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(200,7,'Winter solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(175,4,'Mountains',color='black',fontsize=14,rotation=-17)
        ax0.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax1 = fig.add_subplot(gs[0, 1])
        cf1 = ax1.scatter(azi_scan_fullsky, ele_scan_fullsky, c=data_scan_fullsky[idx_freqs[1],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=8)
        plot_context(ax1)
        ax1.set_xlim(105, 255)
        ax1.set_ylim(-5, 75)
        ax1.text(110, 70, str(int(freqs_MHz[idx_freqs[1]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax1.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False, labelbottom=False)
        ax1.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax2 = fig.add_subplot(gs[1, 0])
        cf2 = ax2.scatter(azi_scan_fullsky, ele_scan_fullsky, c=data_scan_fullsky[idx_freqs[2],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=8)
        ax2.set_xlabel('Azimuth [deg]')
        ax2.set_ylabel('Elevation [deg]')
        plot_context(ax2)
        ax2.set_xlim(105, 255)
        ax2.set_ylim(-5, 75)
        ax2.text(110, 70, str(int(freqs_MHz[idx_freqs[2]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax2.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True)
        ax2.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax2.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        ax3 = fig.add_subplot(gs[1, 1])
        cf3 = ax3.scatter(azi_scan_fullsky, ele_scan_fullsky, c=data_scan_fullsky[idx_freqs[3],:], cmap='gnuplot', vmin= vmin, vmax=vmax, lw=8)
        ax3.set_xlabel('Azimuth [deg]')
        plot_context(ax3)
        ax3.set_xlim(105, 255)
        ax3.set_ylim(-5, 75)
        ax3.text(110, 70, str(int(freqs_MHz[idx_freqs[3]]))+' MHz', fontsize=16, color='black', bbox=dict(facecolor='white', alpha=0.6))
        ax3.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False)
        ax3.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax3.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        cbar = plt.colorbar(cf0, cax=cbar_ax, orientation='horizontal', ticklocation='bottom')
        cbar.set_label('ADU', y=-1.8)  # Move label up by adjusting y

        plt.tight_layout()
        plt.savefig(folder_store+'full-sky-scan_'+str(date_scan_fullsky)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()


        ### Interpolated

        # Create a regular grid to interpolate the data
        x_grid = np.linspace(min(azi_scan_fullsky), max(azi_scan_fullsky), 100)
        y_grid = np.linspace(min(ele_scan_fullsky), max(ele_scan_fullsky), 100)
        X, Y = np.meshgrid(x_grid, y_grid)

        # Interpolate scattered data
        Z1 = griddata((azi_scan_fullsky, ele_scan_fullsky), data_scan_fullsky[idx_freqs[0],:], (X, Y), method='cubic')
        Z2 = griddata((azi_scan_fullsky, ele_scan_fullsky), data_scan_fullsky[idx_freqs[1],:], (X, Y), method='cubic')
        Z3 = griddata((azi_scan_fullsky, ele_scan_fullsky), data_scan_fullsky[idx_freqs[2],:], (X, Y), method='cubic')
        Z4 = griddata((azi_scan_fullsky, ele_scan_fullsky), data_scan_fullsky[idx_freqs[3],:], (X, Y), method='cubic')

        nlevels = 26
        vmin = 138 #ADU_min_TVsat
        vmax = 188 #ADU_max_TVsat


        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(16, 14))

        fig.suptitle('CALLISTO full-sky scanning mode (interpolated) on '+date_scan_fullsky, fontsize=24, y=0.95)

        cbar_ax = fig.add_axes([0.11, 0.06, 0.78, 0.013])

        gs = fig.add_gridspec(2, 2, 
                            top=0.92,     # Higher top margin for title
                            bottom=0.12,  # Higher bottom margin for colorbar
                            left=0.1, 
                            right=0.9,
                            wspace=0.05, 
                            hspace=0.05)

        ax0 = fig.add_subplot(gs[0, 0])
        cf0 = ax0.contourf(X, Y, Z1, cmap='terrain', levels=np.linspace(vmin, vmax, nlevels))
        ax0.set_ylabel('Elevation [deg]')
        plot_context(ax0)
        ax0.set_xlim(105, 255)
        ax0.set_ylim(-5, 75)
        ax0.text(110, 70, str(int(freqs_MHz[idx_freqs[0]]))+' MHz', fontsize=16, color='white', bbox=dict(facecolor='white', alpha=0.3))
        ax0.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelbottom=False)
        ax0.text(222,17,'Geostat. orbit',color='gray',fontsize=14,rotation=-45)
        ax0.text(218,45,'Summer solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(200,7,'Winter solstice',color='black',fontsize=14,rotation=-37)
        ax0.text(175,4,'Mountains',color='black',fontsize=14,rotation=-17)
        ax0.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax1 = fig.add_subplot(gs[0, 1])
        cf1 = ax1.contourf(X, Y, Z2, cmap='terrain', levels=np.linspace(vmin, vmax, nlevels))
        plot_context(ax1)
        ax1.set_xlim(105, 255)
        ax1.set_ylim(-5, 75)
        ax1.text(110, 70, str(int(freqs_MHz[idx_freqs[1]]))+' MHz', fontsize=16, color='white', bbox=dict(facecolor='white', alpha=0.3))
        ax1.tick_params(axis='both', which='major', labelsize=20, direction='in', 
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False, labelbottom=False)
        ax1.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])

        ax2 = fig.add_subplot(gs[1, 0])
        cf2 = ax2.contourf(X, Y, Z3, cmap='terrain', levels=np.linspace(vmin, vmax, nlevels))
        ax2.set_xlabel('Azimuth [deg]')
        ax2.set_ylabel('Elevation [deg]')
        plot_context(ax2)
        ax2.set_xlim(105, 255)
        ax2.set_ylim(-5, 75)
        ax2.text(110, 70, str(int(freqs_MHz[idx_freqs[2]]))+' MHz', fontsize=16, color='white', bbox=dict(facecolor='white', alpha=0.3))
        ax2.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True)
        ax2.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax2.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        ax3 = fig.add_subplot(gs[1, 1])
        cf3 = ax3.contourf(X, Y, Z4, cmap='terrain', levels=np.linspace(vmin, vmax, nlevels))
        ax3.set_xlabel('Azimuth [deg]')
        plot_context(ax3)
        ax3.set_xlim(105, 255)
        ax3.set_ylim(-5, 75)
        ax3.text(110, 70, str(int(freqs_MHz[idx_freqs[3]]))+' MHz', fontsize=16, color='white', bbox=dict(facecolor='white', alpha=0.3))
        ax3.tick_params(axis='both', which='major', labelsize=20, direction='in',
                        length=5, width=2, colors='black', top=True, right=True, labelleft=False)
        ax3.set_xticks([110, 130, 150, 170, 190, 210, 230, 250])
        ax3.set_xticklabels(['110', '130', '150', '170', '190', '210', '230', '250'])

        cbar = plt.colorbar(cf0, cax=cbar_ax, orientation='horizontal', ticklocation='bottom')
        cbar.set_label('ADU', y=-1.8)  # Move label up by adjusting y

        plt.tight_layout()
        plt.savefig(folder_store+'full-sky-scan_interpolated_'+str(date_scan_fullsky)+'.png', dpi=200, 
                    bbox_inches='tight',  # Ensure everything is included in saved figure
                    pad_inches=0.05)      # Add small padding
        plt.close()


        ### If everything runs successfully
        log_execution("Successfully run")
        
    except Exception as e:
        
        ### If an error occurs, log the error message and stack trace
        error_message = f"Error: {str(e)}\n{traceback.format_exc()}"
        log_execution(error_message)
    
if __name__ == '__main__':
    main()
    