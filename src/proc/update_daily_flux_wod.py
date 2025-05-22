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

### Which frequency to use (MHz)
in_freq = 10640

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
    tstop = int((20-Terr)/dT) # start Icold, sec -> pixel
    T2 = int((50-Terr)/dT)
    Icold = np.mean(Dlin[:,tstop:T2],axis=1)  # Icold for each freq, mean over time

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

#************************** MAIN ************************** 
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
    for tstart,tstop in zip(ref_times_start, ref_times_end):
        # concat date and time from ref_times to form proper datetimes
        tstart = datetime.datetime.strptime(this_date.strftime("%Y%m%d ")+tstart,
                                "%Y%m%d %H:%M").replace(tzinfo=datetime.timezone.utc)
        tstop = datetime.datetime.strptime(this_date.strftime("%Y%m%d ")+tstop,
                                "%Y%m%d %H:%M").replace(tzinfo=datetime.timezone.utc)
        df_indexed = df_meteoswiss.loc[(df_meteoswiss["time"] > tstart)
                                & (df_meteoswiss["time"] <= tstop)]

        mean_Thot = df_indexed["temp_degC"].mean() + 273.15
        total_precip = df_indexed["precip_mm"].sum()
        mean_irrad = df_indexed["irrad_W_m2"].mean() 
   
        ### find indices of datetimes_FIT that are in the time interval
        # Get start and end times of fit files
        fit_times_start = files2process[this_date]["times"]
        offset = fit_times_start[-1] - fit_times_start[-2]
        fit_times_stop = fit_times_start[1::]
        # add offset
        fit_times_stop = np.array(list(fit_times_stop) + 
                                [fit_times_stop[-1] + offset])  
        
        isin_timeperiod = np.logical_and(fit_times_start > tstart,
                                         fit_times_stop <= tstop)
        
        if not any(isin_timeperiod): # no file found
            sfu = np.nan
            std = np.nan
            comment = 'TO '
        else:
            indexes = np.where(isin_timeperiod)[0]
            for idx in indexes:
                fitfile = files2process[this_date]['files'][idx]
                sfu = processFIT_returnSFU(fitfile,
                                        this_date, 
                                        mean_Thot, 
                                        total_precip, 
                                        in_freq)
                sfu.extend(sfu)
            comment = get_comment(sfu, total_precip)
            std = integrated_std(sfu)
            sfu = np.mean(sfu)
        
        ### Read the pointing csv file
        YYYY = str(tstart.year)
        filepath_csv_offsets = os.path.join(directory_csv_offsets,
            'pointing-offsets_'+YYYY+'.csv')
        df_offsets = pd.read_csv(filepath_csv_offsets, parse_dates = ["start_scan_UTC"],
                                 date_format  = "%Y-%m-%d %H:%M:%S")    
        start_scan = df_offsets['start_scan_UTC'].dt.tz_localize("UTC")
        
        delta_azi_fit = df_offsets['delta_azi_fit'].values
        delta_ele_fit = df_offsets['delta_ele_fit'].values
        quality_check = df_offsets['quality_check'].values

        # Find index of timesteps that fall in time period
        isin_timeperiod = np.logical_and(start_scan > pd.to_datetime(tstart),
                start_scan <= pd.to_datetime(tstop))

        if not np.any(isin_timeperiod):
            avg_corr = 1
            avg_quality_check = -9999
            
        else:
            # Compute avg off-pointing correction
            dazi = delta_azi_fit[isin_timeperiod]
            dele = delta_ele_fit[isin_timeperiod]
            dist = np.sqrt(dazi**2+dele**2)
            corr = f_offpointing_correction(dist)
            avg_corr = np.mean(corr)
            
            # Compute quality check
            avg_quality_check = np.all(quality_check[isin_timeperiod])

        ### Store the results in the update the csv file in df_flux by adding the new values
        df_flux = pd.read_csv(filepath_csv_flux)
        new_row = pd.DataFrame({'it_start_UT': str(tstart), 
            'it_end_UT': str(tstop), 
            'sfu_10640MHz': sfu, 
            'std_sfu_10640MHz': std, 
            'corr_sfu_10640MHz': sfu/avg_corr,
            'corr_std_sfu_10640MHz': std/avg_corr, 
            'temp_K': mean_Thot, 
            'precip_mm': total_precip, 
            'global_irr_Wm2': mean_irrad,
            'quality_check_pointing': avg_quality_check, 
            'comment': comment}, index=[0])
                    
    #     ### Open the last two CSV files and concatene them
    #     this_date = pd.to_datetime('now')
    #     YYYY = str(this_date.year)
    #     YYYYm1 = str(this_date.year-1)
    #     filename_csv_flux1 = 'daily-flux_%s.csv' % YYYY
    #     filename_csv_flux2 = 'daily-flux_%s.csv' % YYYYm1
    #     filepath_csv_flux1 = os.path.join(directory_output, filename_csv_flux1)
    #     filepath_csv_flux2 = os.path.join(directory_output, filename_csv_flux2)
    #     df_flux1 = pd.read_csv(filepath_csv_flux1)
    #     df_flux2 = pd.read_csv(filepath_csv_flux2)
    #     df_flux = pd.concat([df_flux1, df_flux2])
    #     df_flux = df_flux.sort_values(by='it_start_UT')

    #     ### Divide the df by month and year
    #     now = pd.to_datetime('now')
    #     one_month_ago = now - pd.DateOffset(months=1)
    #     one_year_ago = now - pd.DateOffset(years=1)
    #     df_flux_month = df_flux[pd.to_datetime(df_flux['it_start_UT']) > one_month_ago]
    #     df_flux_year = df_flux[pd.to_datetime(df_flux['it_start_UT']) > one_year_ago]
        
    #     ### Plot the daily fluxes
    #     #plot_daily_fluxes(df_flux_month, directory_output, 'month')
    #     #plot_daily_fluxes(df_flux_year, directory_output, 'year')
        
    #     # If everything runs successfully
    #     log_execution(logfname, "Successfully run")
        
    # except Exception as e:
    #     ### If an error occurs, log the error message and stack trace
    #     error_message = f"Error: {str(e)}\n{traceback.format_exc()}"
    #     log_execution(logfname, error_message)
    