#%%
"""
This script is used to update the daily solar fluxes obtained from CALLISTO. 
By default, it gets the flux at 10.64 GHz. The flux is obtaining by integrating over
an hour between the following time ranges
     - morning:     09:00 - 10:00 UT
     - local noon:  11:00 - 12:00 UT
     - afternoon:   13:00 - 14:00 UT
Afterwards, it also generates a plot of the daily fluxes, over a month and over a year.

Author: Andrea Francesco Battaglia (andrea.francesco.battaglia@irsol.usi.ch)

History:
    - 2025/01/28: [Andrea] Created
    - 2025/02/03: [Andrea] Antenna gain changed from 36.4 to 36.5. All dates have been reprocessed.
    - 2025/02/10: [Andrea] In get_comments, the CI/DM threshold has been lowered from 200 to 170.
                           Not all dates have been reprocessed (only that on 2025/02/09)!!
    - 2025/03/17: [Andrea] The script now can deal with the issue of the absence of calibration measures.
                           This can happen in the case of some tests or for mistakes (see 2025/03/11-2025/03/14)
    - 2025/03/28: [Andrea] We now take into account the correction due to the off-pointing. The new csv files
                           now include a new column, referring to the corrected flux.
                           All dates have been reprocessed (with antenna gain of 36.5). It could be that
                           the fluxes we compared with EOVSA in our first paper are now a bit off.
    - 2025/04/02: [Andrea] Script adapted to allow the calculation of the correction factors with 6 scans rather
                           than only 3, like it was before.
"""

import pandas as pd
import os
import scipy.constants as CON
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import traceback

from datetime import datetime
from astropy.io import fits

import sys
sys.path.append('C:\\xrt\\src\\')

from tools.return_fits2process import *
from tools.utils_fitfile import readfit
from tools.utils import checkdirs

from get_Tcold import get_Tcold
from get_coord_Tcold import get_AziElev_Tcold

from utils import log_execution
import warnings
warnings.filterwarnings('ignore')


#************************** USER INPUTS **************************

### Directory where to store the pointing offsets
directory_csv_offsets = 'C:\\xrt\\output\\solar_scan'

### Frequency to consider for the analysis
freq = 10637 # MHz

### Log filename
logfname = 'C:\\xrt\\output\\data\\proc_daily_flux\\log\\update_daily_flux_log.txt'

### Define periods of time during which to compute sfu
# Data will be averaged within these periods

ref_times_start = ["09:00",
                   "11:00",
                   "13:00"]

ref_times_end = ["10:00",
                 "12:00",
                 "14:00"]

#************************** FUNCTIONS **************************


def return_fits2process_csv(directory_fitfiles, directory_output):
    """
    Description:
        Check which data has not been process yet and return the filenames that need
        to be processed. It returns only dates from the last date processed (in the
        csv files) to the last date available in the raw data. If there is a gap in the
        processed dates, nothing is done!

    Input:
        directory_fitfiles: directory where the raw FITS files are stored
        directory_output: directory where the csv files (with fluxes) are stored

    Output:
        files2process: list of filenames that need to be processed
        dates2process: list of dates of the files that need to be processed
    """

    ### List all csv files in directory_output
    files_csv = [f for f in os.listdir(directory_output) if f.startswith('daily-flux_') and f.endswith('.csv')]
    if len(files_csv) == 0:
        ### In the case no processing has been done yet, we consider from the first FIT file
        ## This is the first date (-1) to consider, according to Marco and Philipp
        latest_date = pd.to_datetime('2024-10-27').date()        
    else:
        ### From the latest csv file, extract the latest date
        files_csv.sort()
        df_latest = pd.read_csv(os.path.join(directory_output, files_csv[-1]))
        latest_date = pd.to_datetime(df_latest['it_start_UT'].iloc[-1]).date()

    ### List all raw FITS files in dir_raw that start with "meteoswiss" and ends with "fit"
    files_raw = [f for f in os.listdir(directory_fitfiles) if f.startswith('meteoswiss') and f.endswith('fit')]
    files_raw.sort()
    ### From files_raw, extract the date of the observation
    dates_raw = [pd.to_datetime(f.split('_')[1]).date() for f in files_raw]
    ### Find the indices of the files in dates_raw that have a date larger than latest_date
    idx = [i for i in range(len(dates_raw)) if dates_raw[i] > latest_date]

    ### Files that have a date larger than latest_date
    files2process = [files_raw[i] for i in idx]
    dates2process = [dates_raw[i] for i in idx]

    ### Remove the elements of the list in dates2process that have the same date
    dates2process = sorted(list(set(dates2process)))
    
    return files2process, dates2process


#**********


def compute_Aeff(gain_dB, freq):

    freq = np.array(freq)
    gain_linear = 10.0**(gain_dB/10.) 
    wavelength = CON.c/freq
    Aeff = wavelength**2 * gain_linear / (4*np.pi)

    return Aeff


#**********


def get_dates_telescope_operations():
    """
    Description:
        From the LogFile of the telescope operations, return the dates of the days
        when the telescope did undergo some operational tests or maintenance.
        These days will be discarded from the analysis.

    Input:
        None

    Output:
        dates: list of dates when the telescope did undergo some operational tests
    """

    path2md = 'C:\\xrt\\LogFile_telescope-operations.md'

    ### Read the markdown file
    with open(path2md, 'r') as file:
        md = file.read()

    ### Split the markdown file into lines
    md_lines = md.split('\n')

    ### Find the indices of the lines with dates
    idx_lines = [i for i, line in enumerate(md_lines) if '202' in line]
    #exclude the first, because it is in the header
    idx_lines = idx_lines[1:]

    ### For each line with a date, find "202" and extract the subsequent date
    dates = []
    for idx in idx_lines:
        date = md_lines[idx].split('202')[1].split(' ')[0]
        dates.append('202'+date)

    return dates


#**********


def processFIT_returnSFU(FIT_file, FIT_date, Thot_K, precip_mm, in_freq):
    """
    Description:
        Process the FITS files and return the daily fluxes in SFU

    Input:
        FIT_file: FITS file to be processed (absolute path)
        FIT_date: date of the FITS file
        Thot_K: temperature of the hot load in K
        precip_mm: precipitation in mm
        in_freq: frequency to be considered in MHz

    Output:
        Ssun: solar flux density observation in SFU
    """

    Terr = 0 # timing error FIT-files in seconds, should be 0 in normal cases
    gain_dB = 36.5 # according to Christian suggestion (see email of 30-12-2024)
    sfu = 1e22

    ### Read the FIT file
    dict_fitfile = readfit(FIT_file)
    data = dict_fitfile['data']         # array of values in raw-date
        # shape of data:   [freqs, time (3600)]
    freqs = dict_fitfile['freqs']       # (original) requency axis in MHz
    timeax = dict_fitfile['timeax']     # time axis in s
    time_start = dict_fitfile['T0']
    dT = dict_fitfile['dT']
    date = dict_fitfile['date'] # e.g. 2024/05/10

    ### Compute Aeff
    Aeff = compute_Aeff(gain_dB, freqs*1e6) 

    ### Convert the raw-dates to dB and then to linear power
    dB = data/255.*2500.0/25.4 # raw-date -> dB
    Dlin = 10.0**(dB/10) # dB -> linear power

    ### Cold calibration
    T1 = int((20-Terr)/dT) # start Icold, sec -> pixel
    T2 = int((50-Terr)/dT)
    Icold = np.mean(Dlin[:,T1:T2],axis=1)  # Icold for each freq, mean over time

    ### Hot calibration
    T3 = int(( 80-Terr)/dT) # start Ihot, sec -> pixel
    T4 = int((110-Terr)/dT)
    Ihot  = np.mean(Dlin[:,T3:T4],axis=1)

    ### Extract the sun measurements
    T5 = int((150-Terr)/dT) # start Isunscan, sec -> pixel
    T6 = int((835-Terr)/dT) # 035: we do not consider the last minute anymore, since there is the solar scan
    sun = Dlin[:,T5:T6]    # for each frequency

    ### Get Tcold for different frequencies (we need to know the elevation)
    time_Tcold = pd.to_datetime(str(FIT_date)+' '+time_start)
    ## if there is no calibration measurements on the day, then get_AziElev_Tcold
    ## we return Ssun = -9999 if this is the case, because without calibration the 
    ## solar flux monitoring does not make any sense in any case.
    try:
        azi_Tcold, elev_Tcold = get_AziElev_Tcold(time_Tcold)
    except:
        msg = 'No calibration found for the file on '+str(FIT_date)+' at '+str(time_start)+'. Solar flux set to -9999'
        log_execution(logfname, msg, plot_separation=False)
        Ssun = [-9999, -9999, -9999, -9999]
        return Ssun
    these_Tcold = get_Tcold(time_Tcold, freqs/1e3, elev_Tcold, precip_mm)

    ### Find the index of the frequency in the input frequency
    idx_freq = np.argmin(np.abs(freqs - in_freq))

    ### Extract the data for the input frequency needed for the solar flux calculation
    this_sun = sun[idx_freq,:]
    #this_red_data = red_data[idx_freq,:]
    this_Icold = Icold[idx_freq]
    this_Ihot = Ihot[idx_freq]
    this_A = Aeff[idx_freq]  # in units of m^2
    #this_wlth = wavelength[idx_freq]
    this_Tcold = these_Tcold[idx_freq]

    ### Get the solar flux in SFU
    Tsun = (Thot_K - this_Tcold) / (this_Ihot - this_Icold) * (this_sun - this_Icold) # + Tcold # Sun temperature in Kelvin
    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
    
    return Ssun


#**********


def get_comment(sfu, precip_mm):
    """
    Description:
        Get the comment for the flux value. It is to know if we potentially have a
        flare or a problem with the calibration.

    Input:
        sfu: flux value in SFU

    Output:
        comment: comment for the flux value
                    - FL/IF: flare or interference
                    - CI/DM: calibration issue or delayed motor
                    - PR: precipitation
    """

    if len(sfu) == 0:
        comment='no data'
    else:    
        comment=''
        if np.max(sfu) > 620: comment = comment+'FL/IF '
        if np.min(sfu) < 170: comment = comment+'CI/DM '
        if precip_mm > 0: comment = comment+'PR '
        if np.mean(sfu) <= -9998 and np.mean(sfu) >= -10000: comment = 'No calibration!'
    
    return comment


#**********


def integrated_std(data):
    """ 
    Description:
        Returns the standard deviation over integrated data
        
    Input:
        data: input data array
        
    Output:
        int_std: integrated standard deviation
    """
    
    ### Integrate over a minute (nominal cadence CALLISTO: 0.25s)
    # Number of points in one minute (60s / 0.25s = 240 points)
    npts = 4*12
    n_times = len(data) // npts
    
    ### Integrate the data over a minute
    int_data = []
    for i in range(n_times):
        start_idx = i * npts
        end_idx = (i + 1) * npts
        mean_int_data = np.mean(data[start_idx:end_idx])
        int_data.append(mean_int_data)
    
    ### Standard deviation of the integrated curve
    int_std = np.std(int_data)
    
    return int_std


#**********


def plot_daily_fluxes(df_flux, directory_output,period):
    """
    Description:
        Plot the daily fluxes over a month and over a year

    Input:
        df_flux: dataframe with the daily fluxes
        directory_output: directory where the plots are stored
        period: 'month' or 'year'

    Output:
        None
    """

    times_all = pd.to_datetime(df_flux['it_start_UT'].values)
    #times_end_all = pd.to_datetime(df_flux['it_end_UT'])
    sfu_all = df_flux['sfu_10640MHz'].values
    std_all = df_flux['std_sfu_10640MHz'].values
    corr_sfu_all = df_flux['corr_sfu_10640MHz'].values
    corr_std_all = df_flux['corr_std_sfu_10640MHz'].values
    precip_all = df_flux['precip_mm'].values
    comment_all = df_flux['comment'].values
    quality_check_all = df_flux['quality_check_pointing'].values
    idx_morning = [i for i in range(len(times_all)) if times_all[i].hour == 9]
    idx_noon = [i for i in range(len(times_all)) if times_all[i].hour == 11]
    idx_afternoon = [i for i in range(len(times_all)) if times_all[i].hour == 13]
    times_morning = times_all[idx_morning]
    times_noon = times_all[idx_noon]
    times_afternoon = times_all[idx_afternoon]
    #times_end_morning = times_end_all[idx_morning]
    #times_end_noon = times_end_all[idx_noon]
    #times_end_afternoon = times_end_all[idx_afternoon]
    sfu_morning = sfu_all[idx_morning]
    sfu_noon = sfu_all[idx_noon]
    sfu_afternoon = sfu_all[idx_afternoon]
    std_morning = std_all[idx_morning]
    std_noon = std_all[idx_noon]
    std_afternoon = std_all[idx_afternoon]
    corr_sfu_morning = corr_sfu_all[idx_morning]
    corr_sfu_noon = corr_sfu_all[idx_noon]
    corr_sfu_afternoon = corr_sfu_all[idx_afternoon]
    corr_std_morning = corr_std_all[idx_morning]
    corr_std_noon = corr_std_all[idx_noon]
    corr_std_afternoon = corr_std_all[idx_afternoon]
    precip_morning = precip_all[idx_morning]
    precip_noon = precip_all[idx_noon]
    precip_afternoon = precip_all[idx_afternoon]
    comment_morning = comment_all[idx_morning]
    comment_noon = comment_all[idx_noon]
    comment_afternoon = comment_all[idx_afternoon]
    quality_check_morning = quality_check_all[idx_morning]
    quality_check_noon = quality_check_all[idx_noon]
    quality_check_afternoon = quality_check_all[idx_afternoon]
    
    ### Get the indices where the pointing quality is bad
    idx_bad_morning = [pp for pp, x in enumerate(quality_check_morning) if x == 0]
    idx_bad_noon = [pp for pp, x in enumerate(quality_check_noon) if x == 0]
    idx_bad_afternoon = [pp for pp, x in enumerate(quality_check_afternoon) if x == 0]
    
    ### Create diagnostics vertical lines
    flares_morning = []
    flares_noon = []
    flares_afternoon = []
    calib_issues_morning = []
    calib_issues_noon = []
    calib_issues_afternoon = []
    precip_morning = []
    precip_noon = []
    precip_afternoon = []
    telescope_operations_morning = []
    telescope_operations_noon = []
    telescope_operations_afternoon = []
    multiple_issues_morning = []
    multiple_issues_noon = []
    multiple_issues_afternoon = []
    for i in range(len(comment_morning)):
        if comment_morning[i] == 'FL/IF ': flares_morning.append(i)
        elif comment_morning[i] == 'CI/DM ': calib_issues_morning.append(i)
        elif comment_morning[i] == 'PR ': precip_morning.append(i)
        elif comment_morning[i] == 'TO ': telescope_operations_morning.append(i)
        elif pd.notna(comment_morning[i]): multiple_issues_morning.append(i)
    for i in range(len(comment_noon)):
        if comment_noon[i] == 'FL/IF ': flares_noon.append(i)
        elif comment_noon[i] == 'CI/DM ': calib_issues_noon.append(i)
        elif comment_noon[i] == 'PR ': precip_noon.append(i)
        elif comment_noon[i] == 'TO ': telescope_operations_noon.append(i)
        elif pd.notna(comment_noon[i]): multiple_issues_noon.append(i)
    for i in range(len(comment_afternoon)):
        if comment_afternoon[i] == 'FL/IF ': flares_afternoon.append(i)
        elif comment_afternoon[i] == 'CI/DM ': calib_issues_afternoon.append(i)
        elif comment_afternoon[i] == 'PR ': precip_afternoon.append(i)
        elif comment_afternoon[i] == 'TO ': telescope_operations_afternoon.append(i)
        elif pd.notna(comment_afternoon[i]): multiple_issues_afternoon.append(i)

    ### Plot the daily fluxes
    
    if period == 'month':
        figsize = (18, 16)
        interval_locator = mdates.DayLocator(interval=3)
    elif period == 'year':
        figsize = (24, 16)
        interval_locator = mdates.MonthLocator(interval=1)
    fig = plt.figure(figsize=figsize)
    plt.rcParams.update({'font.size': 18})
    
    ax1 = fig.add_subplot(311)
    ax1.errorbar(times_morning, corr_sfu_morning, yerr=corr_std_morning, fmt='o', color='gray', label='Morning (09:00 - 10:00 UT), corrected')
    if len(idx_bad_morning) > 0:
        for i in range(len(idx_bad_morning)):
            ax1.scatter(times_morning[idx_bad_morning[i]], corr_sfu_morning[idx_bad_morning[i]], marker='x', color='red', lw=3, s=120)
        ax1.scatter(times_morning[idx_bad_morning[0]], corr_sfu_morning[idx_bad_morning[0]], marker='x', color='red', lw=3, s=120, label='Off-pointing > 0.7 deg')    
    ax1.errorbar(times_morning, sfu_morning, yerr=std_morning, fmt='o', color='black', label='Morning (09:00 - 10:00 UT)')
    #ax1.errorbar(times_noon, sfu_noon, yerr=std_noon, fmt='o', color='black', label='Noon', alpha=0.5)
    #ax1.errorbar(times_afternoon, sfu_afternoon, yerr=std_afternoon, fmt='^', color='black', label='Afternoon', alpha=0.5)
    ax1.set_ylabel('Solar flux density [sfu]')
    ax1.set_title('Daily solar fluxes (morning) over a %s, at 10.64 GHz' % period)
    ax1.xaxis.set_major_locator(interval_locator)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    #ax1.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax1.grid()
    ax1.set_ylim(200, 500)
    ax1.set_xlim(times_all[0]-pd.to_timedelta(10, unit='h'), times_all[-1]+pd.to_timedelta(10, unit='h'))
    for i in range(len(flares_morning)):
        if i == 0:
            ax1.fill_between([times_morning[flares_morning[i]], times_morning[flares_morning[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7, label='Flare or interference')
        else:
            ax1.fill_between([times_morning[flares_morning[i]], times_morning[flares_morning[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7)
    for i in range(len(calib_issues_morning)):
        if i == 0:
            ax1.fill_between([times_morning[calib_issues_morning[i]], times_morning[calib_issues_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='red', alpha=0.7, label='Drop in flux or delayed motor')
        else:
            ax1.fill_between([times_morning[calib_issues_morning[i]], times_morning[calib_issues_morning[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='red', alpha=0.7)
    for i in range(len(precip_morning)):
        if i == 0:    
            ax1.fill_between([times_morning[precip_morning[i]], times_morning[precip_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7, label='Precipitation')
        else:
            ax1.fill_between([times_morning[precip_morning[i]], times_morning[precip_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7)
    for i in range(len(telescope_operations_morning)):
        if i == 0:    
            ax1.fill_between([times_morning[telescope_operations_morning[i]], times_morning[telescope_operations_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7, label='Telescope operation or missing files')
        else:
            ax1.fill_between([times_morning[telescope_operations_morning[i]], times_morning[telescope_operations_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7)
    for i in range(len(multiple_issues_morning)):
        if i == 0:    
            ax1.fill_between([times_morning[multiple_issues_morning[i]], times_morning[multiple_issues_morning[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='black', alpha=0.7, label='Multiple issues', hatch='/')
        else:
            ax1.fill_between([times_morning[multiple_issues_morning[i]], times_morning[multiple_issues_morning[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='black', alpha=0.7, hatch='/')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2 = fig.add_subplot(312)
    #ax2.errorbar(times_morning, sfu_morning, yerr=std_morning, fmt='*', color='black', label='Morning', alpha=0.5)
    ax2.errorbar(times_noon, corr_sfu_noon, yerr=corr_std_noon, fmt='o', color='gray', label='Noon (11:00 - 12:00 UT), corrected')
    if len(idx_bad_noon) > 0:
        for i in range(len(idx_bad_noon)):
            ax2.scatter(times_noon[idx_bad_noon[i]], corr_sfu_noon[idx_bad_noon[i]], marker='x', color='red', lw=3, s=120)
        ax2.scatter(times_noon[idx_bad_noon[0]], corr_sfu_noon[idx_bad_noon[0]], marker='x', color='red', lw=3, s=120, label='Off-pointing > 0.7 deg')
    ax2.errorbar(times_noon, sfu_noon, yerr=std_noon, fmt='o', color='black', label='Noon (11:00 - 12:00 UT)')
    #ax2.errorbar(times_afternoon, sfu_afternoon, yerr=std_afternoon, fmt='^', color='black', label='Afternoon', alpha=0.5)
    ax2.set_ylabel('Solar flux density [sfu]')
    ax2.set_title('Daily solar fluxes (noon) over a %s, at 10.64 GHz' % period)
    ax2.xaxis.set_major_locator(interval_locator)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax2.grid()
    ax2.set_ylim(200, 500)
    ax2.set_xlim(times_all[0]-pd.to_timedelta(10, unit='h'), times_all[-1]+pd.to_timedelta(10, unit='h'))
    for i in range(len(flares_noon)):
        if i == 0:
            ax2.fill_between([times_noon[flares_noon[i]], times_noon[flares_noon[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7, label='Flare or interference')
        else:
            ax2.fill_between([times_noon[flares_noon[i]], times_noon[flares_noon[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7)
    for i in range(len(calib_issues_noon)):
        if i == 0:
            ax2.fill_between([times_noon[calib_issues_noon[i]], times_noon[calib_issues_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='red', alpha=0.7, label='Drop in flux or delayed motor')
        else:
            ax2.fill_between([times_noon[calib_issues_noon[i]], times_noon[calib_issues_noon[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='red', alpha=0.7)
    for i in range(len(precip_noon)):
        if i == 0:    
            ax2.fill_between([times_noon[precip_noon[i]], times_noon[precip_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7, label='Precipitation')
        else:
            ax2.fill_between([times_noon[precip_noon[i]], times_noon[precip_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7)
    for i in range(len(telescope_operations_noon)):
        if i == 0:    
            ax2.fill_between([times_noon[telescope_operations_noon[i]], times_noon[telescope_operations_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7, label='Telescope operation or missing files')
        else:
            ax2.fill_between([times_noon[telescope_operations_noon[i]], times_noon[telescope_operations_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7)
    for i in range(len(multiple_issues_noon)):
        if i == 0:    
            ax2.fill_between([times_noon[multiple_issues_noon[i]], times_noon[multiple_issues_noon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='black', alpha=0.7, label='Multiple issues', hatch='/')
        else:
            ax2.fill_between([times_noon[multiple_issues_noon[i]], times_noon[multiple_issues_noon[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='black', alpha=0.7, hatch='/')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    
    ax3 = fig.add_subplot(313)
    #ax3.errorbar(times_morning, sfu_morning, yerr=std_morning, fmt='*', color='black', label='Morning', alpha=0.5)
    #ax3.errorbar(times_noon, sfu_noon, yerr=std_noon, fmt='o', color='black', label='Noon', alpha=0.5)
    ax3.errorbar(times_afternoon, corr_sfu_afternoon, yerr=corr_std_afternoon, fmt='o', color='gray', label='Afternoon (13:00 - 14:00 UT), corrected')
    if len(idx_bad_afternoon) > 0:
        for i in range(len(idx_bad_afternoon)):
            ax3.scatter(times_afternoon[idx_bad_afternoon[i]], corr_sfu_afternoon[idx_bad_afternoon[i]], marker='x', color='red', lw=3, s=120)
        ax3.scatter(times_afternoon[idx_bad_afternoon[0]], corr_sfu_afternoon[idx_bad_afternoon[0]], marker='x', color='red', lw=3, s=120, label='Off-pointing > 0.7 deg')
    ax3.errorbar(times_afternoon, sfu_afternoon, yerr=std_afternoon, fmt='o', color='black', label='Afternoon (13:00 - 14:00 UT)')
    ax3.set_ylabel('Solar flux density [sfu]')
    ax3.set_title('Daily solar fluxes (afternoon) over a %s, at 10.64 GHz' % period)
    ax3.xaxis.set_major_locator(interval_locator)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    #ax3.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax3.grid()
    ax3.set_ylim(200, 500)
    ax3.set_xlim(times_all[0]-pd.to_timedelta(10, unit='h'), times_all[-1]+pd.to_timedelta(10, unit='h'))
    for i in range(len(flares_afternoon)):
        if i == 0:
            ax3.fill_between([times_afternoon[flares_afternoon[i]], times_afternoon[flares_afternoon[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7, label='Flare or interference')
        else:
            ax3.fill_between([times_afternoon[flares_afternoon[i]], times_afternoon[flares_afternoon[i]]+pd.Timedelta(hours=1)], 
                0, 1000, color='orange', alpha=0.7)
    for i in range(len(calib_issues_afternoon)):
        if i == 0:
            ax3.fill_between([times_afternoon[calib_issues_afternoon[i]], times_afternoon[calib_issues_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='red', alpha=0.7, label='Drop in flux or delayed motor')
        else:
            ax3.fill_between([times_afternoon[calib_issues_afternoon[i]], times_afternoon[calib_issues_afternoon[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='red', alpha=0.7)
    for i in range(len(precip_afternoon)):
        if i == 0:    
            ax3.fill_between([times_afternoon[precip_afternoon[i]], times_afternoon[precip_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7, label='Precipitation')
        else:
            ax3.fill_between([times_afternoon[precip_afternoon[i]], times_afternoon[precip_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='blue', alpha=0.7)
    for i in range(len(telescope_operations_afternoon)):
        if i == 0:    
            ax3.fill_between([times_afternoon[telescope_operations_afternoon[i]], times_afternoon[telescope_operations_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7, label='Telescope operation or missing files')
        else:
            ax3.fill_between([times_afternoon[telescope_operations_afternoon[i]], times_afternoon[telescope_operations_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='green', alpha=0.7)
    for i in range(len(multiple_issues_afternoon)):
        if i == 0:    
            ax3.fill_between([times_afternoon[multiple_issues_afternoon[i]], times_afternoon[multiple_issues_afternoon[i]]+pd.Timedelta(hours=1)],
                0, 1000, color='black', alpha=0.7, label='Multiple issues', hatch='/')
        else:
            ax3.fill_between([times_afternoon[multiple_issues_afternoon[i]], times_afternoon[multiple_issues_afternoon[i]]+pd.Timedelta(hours=1)],
            0, 1000, color='black', alpha=0.7, hatch='/')
    ax3.set_xlabel('Start time: '+str(times_all[0].year)+'-'+str(times_all[0].month)+'-'+str(times_all[0].day))
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(directory_output, 'daily_fluxes_%s.png' % period))
    plt.close()


#**********


def f_offpointing_correction(dist):
    '''
    This is the formula we use to calculate the correction factor for the offpointing.
    In the case of questions, ask Marco Gabella
    
    -----> for the half-power-beam-width we assume 2.8 deg. However, one can get the value
           directly from the csv (hence from the fit; needs to be checked)
    
    dist: distance of the fit from the center of the scan map
    '''
    hpbw = 2.8
    return np.exp(-4*np.log(2)*(dist/hpbw)**2)


#**********


def get_correction_factor(directory_offsets, date):
    '''
    This function calculates and returns the off-pointing correction factor for the flux.
    So far [2024-04-02], we do an average between two consecutive scans to get the final factor.
    This factor needs to be multiplied to the final flux.
    
    
    -----> for the half-power-beam-width we assume 2.8 deg. However, one can get the value
           directly from the csv (hence from the fit; needs to be checked)
    
    '''
    
    date = pd.to_datetime(date)
    
    ### Read the csv file
    YYYY = str(date.year)
    filepath_csv_offsets = os.path.join(directory_offsets, 'pointing-offsets_'+YYYY+'.csv')
    df_offsets = pd.read_csv(filepath_csv_offsets)    
    start_scan = pd.to_datetime(df_offsets['start_scan_UTC'].values)
    delta_azi_fit = df_offsets['delta_azi_fit'].values
    delta_ele_fit = df_offsets['delta_ele_fit'].values
    quality_check = df_offsets['quality_check'].values
    
    ### We need to distinguish the days before and after 2025-04-01, as
    ## from 2025-03-19 to 2025-03-31 only 3 scans were done
    if date > pd.to_datetime('2025-03-31 23:00:00'):
        ### Get the proper times
        time_scan_morning_before = pd.to_datetime(str(date.date()) + ' 08:59:00')
        time_scan_morning_after = pd.to_datetime(str(date.date()) + ' 09:59:00')
        time_scan_noon_before = pd.to_datetime(str(date.date()) + ' 10:59:00')
        time_scan_noon_after = pd.to_datetime(str(date.date()) + ' 11:59:00')
        time_scan_afternoon_before = pd.to_datetime(str(date.date()) + ' 12:59:00')
        time_scan_afternoon_after = pd.to_datetime(str(date.date()) + ' 13:59:00')
        idx_morning_before = np.argmin(np.abs(start_scan-time_scan_morning_before))
        idx_morning_after = np.argmin(np.abs(start_scan-time_scan_morning_after))
        idx_noon_before = np.argmin(np.abs(start_scan-time_scan_noon_before))
        idx_noon_after = np.argmin(np.abs(start_scan-time_scan_noon_after))
        idx_afternoon_before = np.argmin(np.abs(start_scan-time_scan_afternoon_before))
        idx_afternoon_after = np.argmin(np.abs(start_scan-time_scan_afternoon_after))
        
        ### Get the delta azimuth and elevations
        dazi_morning_before = delta_azi_fit[idx_morning_before]
        dele_morning_before = delta_ele_fit[idx_morning_before]
        dazi_morning_after = delta_azi_fit[idx_morning_after]
        dele_morning_after = delta_ele_fit[idx_morning_after]
        dazi_noon_before = delta_azi_fit[idx_noon_before]
        dele_noon_before = delta_ele_fit[idx_noon_before]
        dazi_noon_after = delta_azi_fit[idx_noon_after]
        dele_noon_after = delta_ele_fit[idx_noon_after]
        dazi_afternoon_before = delta_azi_fit[idx_afternoon_before]
        dele_afternoon_before = delta_ele_fit[idx_afternoon_before]
        dazi_afternoon_after = delta_azi_fit[idx_afternoon_after]
        dele_afternoon_after = delta_ele_fit[idx_afternoon_after]
        
        ### Get the quality checks
        quality_check_morning_before = quality_check[idx_morning_before]
        quality_check_morning_after = quality_check[idx_morning_after]
        quality_check_noon_before = quality_check[idx_noon_before]
        quality_check_noon_after = quality_check[idx_noon_after]
        quality_check_afternoon_before = quality_check[idx_afternoon_before]
        quality_check_afternoon_after = quality_check[idx_afternoon_after]
        
        ### Calculate the correction factors (for the formula, ask Marco Gabella)
        dist_morning_before = np.sqrt(dazi_morning_before**2+dele_morning_before**2)
        dist_morning_after = np.sqrt(dazi_morning_after**2+dele_morning_after**2)
        dist_noon_before = np.sqrt(dazi_noon_before**2+dele_noon_before**2)
        dist_noon_after = np.sqrt(dazi_noon_after**2+dele_noon_after**2)
        dist_afternoon_before = np.sqrt(dazi_afternoon_before**2+dele_afternoon_before**2)
        dist_afternoon_after = np.sqrt(dazi_afternoon_after**2+dele_afternoon_after**2)
        c_morning_before = f_offpointing_correction(dist_morning_before)
        c_morning_after = f_offpointing_correction(dist_morning_after)
        c_noon_before = f_offpointing_correction(dist_noon_before)
        c_noon_after = f_offpointing_correction(dist_noon_after)
        c_afternoon_before = f_offpointing_correction(dist_afternoon_before)
        c_afternoon_after = f_offpointing_correction(dist_afternoon_after)
        c_morning = np.mean([c_morning_before, c_morning_after])
        c_noon = np.mean([c_noon_before, c_noon_after])
        c_afternoon = np.mean([c_afternoon_before, c_afternoon_after])
        
        ### quality factor to return
        quality_pointing_morning = 1
        quality_pointing_noon = 1
        quality_pointing_afternoon = 1
        if quality_check_morning_before == 0 or quality_check_morning_after == 0: quality_pointing_morning = 0
        if quality_check_noon_before == 0 or quality_check_noon_after == 0: quality_pointing_noon = 0
        if quality_check_afternoon_before == 0 or quality_check_afternoon_after == 0: quality_pointing_afternoon = 0 
        
    else:
        ### Get the proper times
        time_scan_morning = pd.to_datetime(str(date.date()) + ' 08:59:00')
        time_scan_noon = pd.to_datetime(str(date.date()) + ' 10:59:00')
        time_scan_afternoon = pd.to_datetime(str(date.date()) + ' 12:59:00')
        idx_morning = np.argmin(np.abs(start_scan-time_scan_morning))
        idx_noon = np.argmin(np.abs(start_scan-time_scan_noon))
        idx_afternoon = np.argmin(np.abs(start_scan-time_scan_afternoon))

        ### Get the delta azimuth and elevations
        dazi_morning = delta_azi_fit[idx_morning]
        dele_morning = delta_ele_fit[idx_morning]
        dazi_noon = delta_azi_fit[idx_noon]
        dele_noon = delta_ele_fit[idx_noon]
        dazi_afternoon = delta_azi_fit[idx_afternoon]
        dele_afternoon = delta_ele_fit[idx_afternoon]

        ### Get the quality checks
        quality_pointing_morning = quality_check[idx_morning]
        quality_pointing_noon = quality_check[idx_noon]
        quality_pointing_afternoon = quality_check[idx_afternoon]
        
        ### Calculate the correction factors (for the formula, ask Marco Gabella)
        dist_morning = np.sqrt(dazi_morning**2+dele_morning**2)
        dist_noon = np.sqrt(dazi_noon**2+dele_noon**2)
        dist_afternoon = np.sqrt(dazi_afternoon**2+dele_afternoon**2)
        c_morning = f_offpointing_correction(dist_morning)
        c_noon = f_offpointing_correction(dist_noon)
        c_afternoon = f_offpointing_correction(dist_afternoon)
        
    return c_morning, c_noon, c_afternoon, quality_pointing_morning, quality_pointing_noon, quality_pointing_afternoon


#************************** MAIN **************************


if __name__ == '__main__':

    ### Frequency to be considered in MHz
    in_freq = 10640

    ### Directories
    directory_fitfiles = 'C:\\xrt\\output\\data\\raw\\FITfiles'
    directory_output = 'C:\\xrt\\output\\data\\proc_daily_flux'
    directory_meteoswiss = 'C:\\xrt\\output\\data\\meteoswiss'
    directory_offsets = 'C:\\xrt\\output\\solar_scan'
    
    checkdirs([directory_fitfiles, directory_output])

    ### This try/except is used to track errors and prinf them in the log file
    try:
            
        ### Get the filenames and dates of the files to process
        filenames2process, dates2process = return_fits2process_csv(directory_fitfiles, directory_output)

        ### Get the dates of the telescope operations
        dates_telescope_operations = get_dates_telescope_operations()

        ### Loop on all the days to process
        len_dates = len(dates2process)
        if len(dates2process) > 0:
            for this_date in dates2process:

                ### Find the meteoswiss CSV file and open it
                YYYY = str(this_date.year)
                MM = str(this_date.month) if this_date.month >= 10 else '0'+str(this_date.month)
                DD = str(this_date.day) if this_date.day >= 10 else '0'+str(this_date.day)
                this_date_csv = YYYY+'-'+MM+'-'+DD
                filename_csv = 'meteoswiss_%s.csv' % this_date_csv 
                filepath_csv = os.path.join(directory_meteoswiss, filename_csv)
                df_meteoswiss = pd.read_csv(filepath_csv, sep=',')
                time_meteoswiss = pd.to_datetime(df_meteoswiss['time'].values)
                Thot_meteoswiss = df_meteoswiss['temp_degC'].values
                precip_meteoswiss = df_meteoswiss['precip_mm'].values
                irrad_meteoswiss = df_meteoswiss['irrad_W_m2'].values

                ### Find the csv with the daily fluxes
                filename_csv_flux = 'daily-flux_%s.csv' % YYYY
                filepath_csv_flux = os.path.join(directory_output, filename_csv_flux)
                #print(filepath_csv_flux)
                # if it does not exist, create it
                if os.path.exists(filepath_csv_flux):
                    df_flux = pd.read_csv(filepath_csv_flux)
                else:
                    df_flux = pd.DataFrame(columns=['it_start_UT', 'it_end_UT', 'sfu_10640MHz', 'std_sfu_10640MHz', 'corr_sfu_10640MHz', 'corr_std_sfu_10640MHz', 'temp_K', 'precip_mm', 'global_irr_Wm2', 'quality_check_pointing', 'comment'])
                    df_flux.to_csv(filepath_csv_flux, index=False)

                ### Find the indices of the measures in the morning, noon and afternoon
                idx_morning = [i for i in range(len(time_meteoswiss)) if time_meteoswiss[i].hour == 9]
                idx_noon = [i for i in range(len(time_meteoswiss)) if time_meteoswiss[i].hour == 11]
                idx_afternoon = [i for i in range(len(time_meteoswiss)) if time_meteoswiss[i].hour == 13]

                ### Mean weather values for morning, noon and afternoon
                mean_Thot_morning = Thot_meteoswiss[idx_morning].mean() + 273.15
                mean_Thot_noon = Thot_meteoswiss[idx_noon].mean() + 273.15
                mean_Thot_afternoon = Thot_meteoswiss[idx_afternoon].mean() + 273.15
                total_precip_morning = precip_meteoswiss[idx_morning].sum()
                total_precip_noon = precip_meteoswiss[idx_noon].sum()
                total_precip_afternoon = precip_meteoswiss[idx_afternoon].sum()
                mean_irrad_morning = irrad_meteoswiss[idx_morning].mean()
                mean_irrad_noon = irrad_meteoswiss[idx_noon].mean()
                mean_irrad_afternoon = irrad_meteoswiss[idx_afternoon].mean()

                ### Get all the filenames of the FIT files for this date
                date_FIT = str(this_date).replace('-', '')
                filenames_this_date = [f for f in filenames2process if f.split('_')[1] == date_FIT]

                ### Datetimes of the FIT files
                datetimes_FIT = [pd.to_datetime(f.split('_')[1]+' '+f.split('_')[2]) for f in filenames_this_date]

                ### find indices of datetimes_FIT that are in the morning, noon and afternoon
                idx_morning_FIT = [i for i in range(len(datetimes_FIT)) if datetimes_FIT[i].hour == 9]
                idx_noon_FIT = [i for i in range(len(datetimes_FIT)) if datetimes_FIT[i].hour == 11]
                idx_afternoon_FIT = [i for i in range(len(datetimes_FIT)) if datetimes_FIT[i].hour == 13]

                ### Loop on all files at the different times of the day
                # morning
                if len(idx_morning_FIT) != 4 or this_date in dates_telescope_operations:
                    print()
                    print('Warning: telescope operation or missing FIT files on %s (morning)' % this_date)
                    print('NaNs will be stored in the output file')
                    sfu_morning = np.nan
                    std_morning = np.nan
                    comment_morning = 'TO '
                else:
                    sfu_morning = []
                    for i in idx_morning_FIT:
                        path_fitfile = os.path.join(directory_fitfiles, filenames_this_date[i])
                        sfu_m = processFIT_returnSFU(path_fitfile, this_date, mean_Thot_morning, total_precip_morning, in_freq)
                        sfu_morning.extend(sfu_m)
                    comment_morning = get_comment(sfu_morning, total_precip_morning)
                    std_morning = integrated_std(sfu_morning)
                    sfu_morning = np.mean(sfu_morning)

                # noon
                if len(idx_noon_FIT) != 4 or this_date in dates_telescope_operations:
                    print()
                    print('Warning: telescope operation or missing FIT files on %s (noon)' % this_date)
                    print('NaNs will be stored in the output file')
                    sfu_noon = np.nan
                    std_noon = np.nan
                    comment_noon = 'TO '
                else:        
                    sfu_noon = []
                    for i in idx_noon_FIT:
                        path_fitfile = os.path.join(directory_fitfiles, filenames_this_date[i])
                        sfu_n = processFIT_returnSFU(path_fitfile, this_date, mean_Thot_noon, total_precip_noon, in_freq)
                        sfu_noon.extend(sfu_n)
                    comment_noon = get_comment(sfu_noon, total_precip_noon)
                    std_noon = integrated_std(sfu_noon)
                    sfu_noon = np.mean(sfu_noon)

                # afternoon
                if len(idx_afternoon_FIT) != 4 or this_date in dates_telescope_operations:
                    print()
                    print('Warning: telescope operation or missing FIT files on %s (afternoon)' % this_date)
                    print('NaNs will be stored in the output file')
                    sfu_afternoon = np.nan
                    std_afternoon = np.nan
                    comment_afternoon = 'TO '
                else:
                    sfu_afternoon = []
                    for i in idx_afternoon_FIT:
                        path_fitfile = os.path.join(directory_fitfiles, filenames_this_date[i])
                        sfu_a = processFIT_returnSFU(path_fitfile, this_date, mean_Thot_afternoon, total_precip_afternoon, in_freq)
                        sfu_afternoon.extend(sfu_a)
                    comment_afternoon = get_comment(sfu_afternoon, total_precip_afternoon)
                    std_afternoon = integrated_std(sfu_afternoon)
                    sfu_afternoon = np.mean(sfu_afternoon)
                
                
            
#%%
                ### Get the correction factors
                # the solar scans started only on 2025-03-19
                if pd.to_datetime(this_date) > pd.to_datetime('2025-03-18 18:00:00'):
                    c_morning, c_noon, c_afternoon, quality_pointing_morning, quality_pointing_noon, quality_pointing_afternoon = get_correction_factor(directory_offsets, this_date)
                else:
                    c_morning = 1
                    c_noon = 1
                    c_afternoon = 1
                    quality_pointing_morning = 1
                    quality_pointing_noon = 1
                    quality_pointing_afternoon = 1
                
                ### Store the results in the update the csv file in df_flux by adding the new values
                df_flux = pd.read_csv(filepath_csv_flux)
                new_row = pd.DataFrame({'it_start_UT': str(this_date)+' 09:00:00', 
                                        'it_end_UT': str(this_date)+' 10:00:00', 
                                        'sfu_10640MHz': sfu_morning, 
                                        'std_sfu_10640MHz': std_morning, 
                                        'corr_sfu_10640MHz': sfu_morning/c_morning,
                                        'corr_std_sfu_10640MHz': std_morning/c_morning, 
                                        'temp_K': mean_Thot_morning, 
                                        'precip_mm': total_precip_morning, 
                                        'global_irr_Wm2': mean_irrad_morning,
                                        'quality_check_pointing': quality_pointing_morning, 
                                        'comment': comment_morning}, index=[0])
                df_flux = pd.concat([df_flux,new_row], ignore_index=True)
                new_row = pd.DataFrame({'it_start_UT': str(this_date)+' 11:00:00',
                                            'it_end_UT': str(this_date)+' 12:00:00',
                                            'sfu_10640MHz': sfu_noon,
                                            'std_sfu_10640MHz': std_noon,
                                            'corr_sfu_10640MHz': sfu_noon/c_noon,
                                            'corr_std_sfu_10640MHz': std_noon/c_noon,
                                            'temp_K': mean_Thot_noon,
                                            'precip_mm': total_precip_noon,
                                            'global_irr_Wm2': mean_irrad_noon,
                                            'quality_check_pointing': quality_pointing_noon,
                                            'comment': comment_noon}, index=[0])
                df_flux = pd.concat([df_flux,new_row], ignore_index=True)
                new_row = pd.DataFrame({'it_start_UT': str(this_date)+' 13:00:00',
                                            'it_end_UT': str(this_date)+' 14:00:00',
                                            'sfu_10640MHz': sfu_afternoon,
                                            'std_sfu_10640MHz': std_afternoon,
                                            'corr_sfu_10640MHz': sfu_afternoon/c_afternoon,
                                            'corr_std_sfu_10640MHz': std_afternoon/c_afternoon,
                                            'temp_K': mean_Thot_afternoon,
                                            'precip_mm': total_precip_afternoon,
                                            'global_irr_Wm2': mean_irrad_afternoon,
                                            'quality_check_pointing': quality_pointing_afternoon,
                                            'comment': comment_afternoon}, index=[0])
                df_flux = pd.concat([df_flux,new_row], ignore_index=True)
                df_flux = df_flux.sort_values(by='it_start_UT')
                #df_flux.to_csv(filepath_csv_flux, index=False)
        
        ### Open the last two CSV files and concatene them
        this_date = pd.to_datetime('now')
        YYYY = str(this_date.year)
        YYYYm1 = str(this_date.year-1)
        filename_csv_flux1 = 'daily-flux_%s.csv' % YYYY
        filename_csv_flux2 = 'daily-flux_%s.csv' % YYYYm1
        filepath_csv_flux1 = os.path.join(directory_output, filename_csv_flux1)
        filepath_csv_flux2 = os.path.join(directory_output, filename_csv_flux2)
        df_flux1 = pd.read_csv(filepath_csv_flux1)
        df_flux2 = pd.read_csv(filepath_csv_flux2)
        df_flux = pd.concat([df_flux1, df_flux2])
        df_flux = df_flux.sort_values(by='it_start_UT')

        ### Divide the df by month and year
        now = pd.to_datetime('now')
        one_month_ago = now - pd.DateOffset(months=1)
        one_year_ago = now - pd.DateOffset(years=1)
        df_flux_month = df_flux[pd.to_datetime(df_flux['it_start_UT']) > one_month_ago]
        df_flux_year = df_flux[pd.to_datetime(df_flux['it_start_UT']) > one_year_ago]
        
        ### Plot the daily fluxes
        #plot_daily_fluxes(df_flux_month, directory_output, 'month')
        #plot_daily_fluxes(df_flux_year, directory_output, 'year')
        
        # If everything runs successfully
        log_execution(logfname, "Successfully run")
        
    except Exception as e:
        ### If an error occurs, log the error message and stack trace
        error_message = f"Error: {str(e)}\n{traceback.format_exc()}"
        log_execution(logfname, error_message)
    