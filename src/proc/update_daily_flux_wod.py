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
import glob
import scipy.constants as CON
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import traceback
import datetime
from astropy.io import fits

import sys
sys.path.append('C:\\xrt\\src\\')

from tools.return_fits2process import *
from tools.utils_fitfile import readfit
from tools.utils import checkdirs
from tools.get_Tcold import get_Tcold
from tools.get_coord_Tcold import get_AziElev_Tcold

from utils import log_execution
from utils import get_fit_files
from utils import get_meteoswiss_data

import warnings
warnings.filterwarnings('ignore')


#************************** USER INPUTS **************************

### Directory where to store the pointing offsets
directory_csv_offsets = 'C:\\xrt\\output\\solar_scan'

### Directory where the dailys sfu values are saved
directory_output = 'C:\\xrt\\output\\data\\proc_daily_flux'

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


def return_fits2process_csv(directory_output):
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
    files_csv = glob.glob(os.path.join(directory_output, "daily-flux*.csv"))
    
    latest_date = datetime.datetime(2024, 10, 27) # default if no file available
    if len(files_csv):
        ### In the case no processing has been done yet, we consider from the first FIT file
        ## This is the first date (-1) to consider, according to Marco and Philipp
        df_latest = pd.read_csv(files_csv[-1])
        latest_date = pd.to_datetime(df_latest['it_start_UT'].iloc[-1])

    ctime = latest_date
    # loop on all days from latest_date to now
    files2process = {}
    while ctime < datetime.datetime.now():
        temp = get_fit_files(ctime.strftime("%Y%m%d"))
        # key is first date of day, value is dict with keys
        # times (all times) and files (all filenames)
        files2process[temp[1][0]] = {"files": temp[0], "times": temp[1]}
        ctime += datetime.timedelta(days = 1)
    
    return files2process


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

### This try/except is used to track errors and prinf them in the log file
    
### Get the filenames and dates of the files to process
files2process = return_fits2process_csv(directory_output)

### Get the dates of the telescope operations
dates_telescope_operations = get_dates_telescope_operations()


if not len(files2process):
    sys.exit()
    
### Loop on all the days to process
len_dates = len(files2process)

for this_date in files2process:
    df_meteoswiss = get_meteoswiss_data(this_date)

    ### Find the csv with the daily fluxes
    filename_csv_flux = f'daily-flux_{this_date.strftime("%Y")}.csv'
    filepath_csv_flux = os.path.join(directory_output, filename_csv_flux)
    
    # if it does not exist, create it
    if os.path.exists(filepath_csv_flux):
        df_flux = pd.read_csv(filepath_csv_flux)
    else:
        df_flux = pd.DataFrame(columns=['it_start_UT', 'it_end_UT', 'sfu_10640MHz', 'std_sfu_10640MHz', 'corr_sfu_10640MHz', 'corr_std_sfu_10640MHz', 'temp_K', 'precip_mm', 'global_irr_Wm2', 'quality_check_pointing', 'comment'])
        df_flux.to_csv(filepath_csv_flux, index=False)

    # Loop on all time periods
    for t0,t1 in zip(ref_times_start, ref_times_end):
        # concat date and time from ref_times to form proper datetimes
        t0 = datetime.datetime.strptime(this_date.strftime("%Y%m%d ")+t0,
                                        "%Y%m%d %H:%M")
        t1 = datetime.datetime.strptime(this_date.strftime("%Y%m%d ")+t1,
                                        "%Y%m%d %H:%M")
        df_indexed = df_meteoswiss.loc[(df_meteoswiss["time"] > t0)
                                & (df_meteoswiss["time"] <= t1)]

        mean_Thot = df_indexed["temp_degC"].mean() + 273.15
        total_precip = df_indexed["precip_mm"].sum()
        mean_irrad = df_indexed["irrad_W_m2"].mean() 
   
        ### find indices of datetimes_FIT that are in the time interval
        fit_times_day = files2process[this_date]["times"]
        idx = np.logical_and(fit_times_day > t0, fit_times_day <= t1)
        
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
    