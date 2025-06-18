# -*- coding: utf-8 -*-

description = """
Description:

    Read FIT-files from local folder. This is an extension of Philipp's script 
    "read_fitfiles.py" but with the addition of processing all frequencies.

    @author: Andrea F. Battaglia

    History:
        - 2024/12/20 [Andrea F. Battaglia]: first version
        - 2024/12/30 [Andrea F. Battaglia]: due to the cross-calibration analysis of
            CALLISTO with EOVSA, we changed the antenna gain_dB from 36 to 36.4
        - 2025/02/03 [Andrea F. Battaglia]: Antenna gain changed from 36.4 to 36.3. 
            Pay attention, because: 
                !!!!!!!! NO DATES HAVE BEEN REPROCESSED !!!!!!!!!!!!!!!!!!!!!
        - 2025/03/15 [Andrea F. Battaglia]: After we started getting night-scan data, 
            the script crashed. The problem has been solved.
"""

import numpy as np
import scipy.constants as CON
import argparse
import glob
import os
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import datetime

from tools.return_fits2process import *
from tools.utils_fitfile import readfit
from tools.utils import checkdirs

import sys
sys.path.append('C:\\xrt\\src\\tools')

from get_Tcold import get_Tcold
from get_coord_Tcold import get_AziElev_Tcold

import warnings
warnings.filterwarnings('ignore')

def compute_Aeff(gain_dB, freq):

    freq = np.array(freq)
    gain_linear = 10.0**(gain_dB/10.) 
    wavelength = CON.c/freq
    Aeff = wavelength**2 * gain_linear / (4*np.pi)

    return Aeff

def main():
    
    # Adding argparser allows to add arguments when running the code from the
    # terminal. Below, two arguments are added:
    #   - to allow the printing of all processing steps
    #   - for the selection of the day
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=description, epilog='')
    
    #parser.add_argument(
    #    '-verbose', action='store_true', required=False, default=False,
    #    help='Print more information like results (default: %(default))')
    
    parser.add_argument(
        '-day', type=str, metavar=str.__name__, 
        action='store', required=False, default=None,
        help='Choose day to analyze (default: analyze all days)')
    
    
    
    args = parser.parse_args()
    #verbose = args.verbose
    day = args.day
    #Thot = args.Thot

    # --------------------------
    # Definitions (to be revised and adjusted!)
    # --------------------------

    # Tcold = 12. # from literature
    Terr = 0 # timing error FIT-files in seconds, should be 0 in normal cases
    #if not Thot:
    #    Thot = 273.15 + 5 # 11.0 # hard coded - it should be taken from the outdoor sensor
    
    #gain_dB = 36 # according to data sheet of TRIAX TDS 65 A
    gain_dB = 36.5 # according to Christian suggestion (see email of 30-12-2024)
    sfu = 1e22

    #freq_MHz = 11075 # in case of 0.25 sec FIT-file
    #dfreq_MHz = 50 # MHz
    # [0081]=11056.562,0
    # [0082]=11064.188,0
    # [0083]=11071.812,0
    # [0084]=11079.438,0
    # [0085]=11087.062,0
    # [0086]=11094.688,0
    
    # --------------------------
    # Define files to be read
    # Create output directory 
    # --------------------------

    #filename_output = 'solarflux.csv'
    filename_output_pkl = 'processed-obs'
    directory_fitfiles = 'C:\\xrt\\output\\data\\raw\\FITfiles'
    directory_output = 'C:\\xrt\\output\\data\\proc'
    directory_meteoswiss = 'C:\\xrt\\output\\data\\meteoswiss'
    checkdirs([directory_fitfiles, directory_output])

    ## If the user run the script with the argument "-day" activated, then enter this
    if day:
        if day == 'all':
            ## If "-day all" is run, then reprocess all files
            filename_format = 'meteoswiss_*_*_01.fit'
            filepath_format = os.path.join(directory_fitfiles, filename_format)
            filepaths = sorted(glob.glob(filepath_format))
        elif day == 'missing':
            ## Only process the files not processed yet (missing processed files)
            ## Detect which files have not been processed yet
            files2process = return_fits2process(directory_fitfiles, directory_output)
            filepaths = [os.path.join(directory_fitfiles, f) for f in files2process]
            filepaths.sort()
        else:
            ## else, only the selected day
            filename_format = 'meteoswiss_%s_*_01.fit' % day #MyFile = 'SWISS-METEO_20231127_*_01.fit'
            filepath_format = os.path.join(directory_fitfiles, filename_format)
            filepaths = sorted(glob.glob(filepath_format))

            ## This is to avoid getting the night-scan data, otherwise the script crashes
            datetimes_files = [pd.to_datetime(os.path.basename(file).split('_')[1] + ' ' + os.path.basename(file).split('_')[2]) for file in filepaths]
            idx_mintime = np.argmin(np.abs(pd.to_datetime(day + ' 08:30') - np.array(datetimes_files)))
            idx_maxtime = np.argmin(np.abs(pd.to_datetime(day + ' 14:30') - np.array(datetimes_files)))
            filepaths = filepaths[idx_mintime:idx_maxtime+1]
        
    else:
        ### If day is not called, then only process today's files
        today = datetime.date.today()
        today_str = today.strftime("%Y%m%d")
        filename_format = 'meteoswiss_%s_*_01.fit' % today_str
        filepath_format = os.path.join(directory_fitfiles, filename_format)
        filepaths = sorted(glob.glob(filepath_format))
        
        ## This is to avoid getting the night-scan data, otherwise the script crashes
        datetimes_files = [pd.to_datetime(os.path.basename(file).split('_')[1] + ' ' + os.path.basename(file).split('_')[2]) for file in filepaths]
        idx_mintime = np.argmin(np.abs(pd.to_datetime(today_str + ' 08:30') - np.array(datetimes_files)))
        idx_maxtime = np.argmin(np.abs(pd.to_datetime(today_str + ' 14:30') - np.array(datetimes_files)))
        filepaths = filepaths[idx_mintime:idx_maxtime+1]

    if len(filepaths) < 1:
        print('All raw data has been processed (or no FIT-file has been found). Return without processing')
        return

    ## Create the folder where to store the processed data
    dates = sorted(list(set([os.path.basename(path)[11:19] for path in filepaths])))
    for date in dates:
        directory_output_day = os.path.join(directory_output, date)
        if not os.path.isdir(directory_output_day):
            print('Create directory: %s' % directory_output_day)
            os.makedirs(directory_output_day)
        #print('Create file: %s' % os.path.join(directory_output_day, filename_output))

    # --------------------------
    # Process files
    # --------------------------

    ## We loop on all different dates
    ## Then we loop on all files on the same dates
    for this_date in dates:
        
        
        # find the meteoswiss CSV file and open it
        this_date_csv = this_date[0:4]+'-'+this_date[4:6]+'-'+this_date[6:8]
        filename_csv = 'meteoswiss_%s.csv' % this_date_csv 
        filepath_csv = os.path.join(directory_meteoswiss, filename_csv)
        df_meteoswiss = pd.read_csv(filepath_csv, sep=',')
        time_meteoswiss = pd.to_datetime(df_meteoswiss['time'].values).tz_localize(None)
        Thot_meteoswiss = df_meteoswiss['temp_degC'].values
        precip_meteoswiss = df_meteoswiss['precip_mm'].values
        #irrad_meteoswiss = df_meteoswiss['irrad_W_m2'].values
        
        
        # if we did some telescope operation on a particular day, the code is likely
        # to fail. This is the reason why we do this, to not stop the iteration
        try:


            ## Get the indices of the files corresponding to this date
            idx_files = [i for i, a in enumerate(filepaths) if '_'+this_date in a]

            ## Open the first file to get the size of the arrays
            #dict_fitfile = readfit(filepath)
            #data = dict_fitfile['data']         # array of values in raw-date
            #freqs = dict_fitfile['freqs']       # (original) requency axis in GHz

            ## We initialize some variables that then will be stored in a pkl file.
            # These variables are useful to diagnose the observations and assess the calibration
            # status during a particular day.
            raw_data           = []     # raw data in linear power
            reduced_data       = []     # reduced data in sfu (including times of Thot and Tcold)
            solar_obs          = []     # reduced solar data (background subtracted) in sfu
            bkg_obs            = []     # background in sfu
            freq_axis_MHz      = []     # frequency axis in MHz
            start_times        = []     # start time of each file
            #idx_Thot           = []     # indices to get the Thot measurements
            #idx_Tcold          = []     # indices to get the Tcold measurements
            #lc_sun_linpow      = []
            sun_linpow         = []     # Solar observations in linear power
            Thot_linpow        = []     # Thot in linear power
            Tcold_linpow       = []     # Tcold in linear power
            time_axis_raw_s    = []     # time axis raw_data
            time_axis_sun_s    = []     # time axis of solar_obs (solar observations)
            time_axis_Thot_s   = []     # time axis Thot
            time_axis_Tcold_s  = []     # time axis Tcold
            Tcold_freq         = []     # Tcold for each frequency

            dt_time_axis = 0
            ## now loop on the files
            for idx_file in idx_files:
                
                this_file = filepaths[idx_file]

                ## Read the fit files
                dict_fitfile = readfit(this_file)
                data = dict_fitfile['data']         # array of values in raw-date
                    # shape of data:   [freqs, time (3600)]
                freqs = dict_fitfile['freqs']       # (original) requency axis in MHz
                timeax = dict_fitfile['timeax']     # time axis in s
                time_start = dict_fitfile['T0']
                dT = dict_fitfile['dT']
                date = dict_fitfile['date'] # e.g. 2024/05/10
                date_stringformat = date[0:4] + date[5:7] + date[8:] # 20240510
                directory_output_day = os.path.join(directory_output, date_stringformat)
                nfreqs = len(freqs)

                ## Increase the time axis
                new_time_axis = timeax + dt_time_axis

                ## Compute Aeff
                Aeff = compute_Aeff(gain_dB, freqs*1e6) 
                gain_linear = 10.0**(gain_dB/10.) 
                wavelength = CON.c/(np.array(freqs)*1e6)
                #Aeff = wavelength**2 * gain_linear / (4*np.pi)

                ## Convert the raw-dates to dB and then to linear power
                dB = data/255.*2500.0/25.4 # raw-date -> dB
                Dlin = 10.0**(dB/10) # dB -> linear power
                
                ## Average the linear power array over a freqency band
                #LC = np.mean(Dlin[freqidx_min:freqidx_max,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files with 100 channels
                
                ## Cold calibration
                T1 = int((20-Terr)/dT) # start Icold, sec -> pixel
                T2 = int((50-Terr)/dT)
                Icold = np.mean(Dlin[:,T1:T2],axis=1)  # Icold for each freq, mean over time
                Icold_long = Dlin[:,T1:T2]
                these_times_Tcold = new_time_axis[T1:T2]

                ## Hot calibration
                T3 = int(( 80-Terr)/dT) # start Ihot, sec -> pixel
                T4 = int((110-Terr)/dT)
                Ihot  = np.mean(Dlin[:,T3:T4],axis=1)
                Ihot_long = Dlin[:,T3:T4]
                these_times_Thot = new_time_axis[T3:T4]

                ## Extract the sun measurements
                #T5 = int((170-Terr)/dT) # start Isunscan, sec -> pixel
                T5 = int((150-Terr)/dT) # start Isunscan, sec -> pixel
                T6 = int((900-Terr)/dT) # 900
                red_data = np.zeros((nfreqs,len(timeax)))
                sun = Dlin[:,T5:T6]    # for each frequency
                #print(sun.shape)
                #print(sun_data.shape)
                #print(len(timeax))
                red_data[:,T5:T6] = sun
                these_times_sun = new_time_axis[T5:T6]
                nsun_time = len(red_data[0,:])

                ## Get Tcold for different frequencies (we need to know the elevation)
                time_Tcold = pd.to_datetime(this_date+' '+time_start)
                #print(time_Tcold)
                idx_precip = np.argmin(np.abs(time_meteoswiss - time_Tcold))
                azi_Tcold, elev_Tcold = get_AziElev_Tcold(time_Tcold)
                these_Tcold = get_Tcold(time_Tcold, freqs/1e3, elev_Tcold,precip_meteoswiss[idx_precip])
                
                ## Get Thot from the measurements of meteoswiss
                diff_time = [abs(time_Tcold - ttt) for ttt in time_meteoswiss]
                idx_Thot = np.argmin(diff_time)
                Thot = Thot_meteoswiss[idx_Thot] + 273.15  # K
                #print(Thot)
                
                #*******************************************************#
                # The following loop can be implemented more efficiently
                # by doing calculation between matrices and not in loop
                #*******************************************************#
                tmp_solar_obs = np.zeros((nfreqs, len(these_times_sun)))
                tmp_red_data = np.zeros((nfreqs, nsun_time))
                tmp_bkg_obs = np.zeros((nfreqs, nsun_time))
                ## Loop on all frequencies
                for idx_freq in range(nfreqs):
                    
                    this_sun = sun[idx_freq,:]
                    this_red_data = red_data[idx_freq,:]
                    this_Icold = Icold[idx_freq]
                    this_Ihot = Ihot[idx_freq]
                    this_A = Aeff[idx_freq]  # in units of m^2
                    this_wlth = wavelength[idx_freq]
                    this_Tcold = these_Tcold[idx_freq]
                    
                    '''
                    ## Philipp's version
                    MyFluxLimit = 0.97 # take only values above 97% of peak flux
                    if len(this_sun) > 0:
                        Isun_max = np.max(this_sun)
                        Isun_mean_above_max = np.mean(this_sun[this_sun>MyFluxLimit*Isun_max])
                    else:
                        Isun_max = np.nan
                        Isun_mean_above_max = np.nan
                    
                    Tsun = (Thot - Tcold) / (this_Ihot - this_Icold) * (Isun_mean_above_max - this_Icold) + Tcold # Antenna temperature in Kelvin
                    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
                    Sback = (2*CON.Boltzmann*Tcold)/(this_wlth**2*gain_linear)*4*np.pi * sfu
                    Ssun_minus_Sback = Ssun - Sback
                    tmp_bkg_obs[idx_freq, :] = Sback
                    tmp_solar_obs[idx_freq, :] = Ssun_minus_Sback
                    '''

                    ### Old formulas with Tcold in Tsun
                    '''
                    Tsun = (Thot - Tcold) / (this_Ihot - this_Icold) * (this_sun - this_Icold) + Tcold # Antenna temperature in Kelvin
                    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
                    Sback = (2*CON.Boltzmann*Tcold)/(this_wlth**2*gain_linear)*4*np.pi * sfu
                    Ssun_minus_Sback = Ssun - Sback
                    tmp_bkg_obs[idx_freq, :] = Sback
                    tmp_solar_obs[idx_freq, :] = Ssun_minus_Sback

                    Tsun = (Thot - Tcold) / (this_Ihot - this_Icold) * (this_red_data - this_Icold) + Tcold # Antenna temperature in Kelvin
                    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
                    Ssun_minus_Sback = Ssun - Sback
                    tmp_red_data[idx_freq, :] = Ssun_minus_Sback
                    '''
                    
                    ### New way to test, without Tcold in the formula
                    Tsun = (Thot - this_Tcold) / (this_Ihot - this_Icold) * (this_sun - this_Icold) # + Tcold # Sun temperature in Kelvin
                    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
                    Sback = (2*CON.Boltzmann*this_Tcold)/(this_wlth**2*gain_linear)*4*np.pi * sfu
                    #Ssun_minus_Sback = Ssun - Sback
                    tmp_bkg_obs[idx_freq, :] = Sback
                    tmp_solar_obs[idx_freq, :] = Ssun

                    Tsun = (Thot - this_Tcold) / (this_Ihot - this_Icold) * (this_red_data - this_Icold) # + Tcold # Antenna temperature in Kelvin
                    Ssun = 2*CON.Boltzmann*Tsun/this_A*sfu
                    #Ssun_minus_Sback = Ssun - Sback
                    tmp_red_data[idx_freq, :] = Ssun # _minus_Sback
                    
                ## Store the variables in loop
                if idx_file == idx_files[0]:
                    raw_data = Dlin
                    solar_obs = tmp_solar_obs
                    reduced_data = tmp_red_data
                    bkg_obs = tmp_bkg_obs
                    sun_linpow = sun
                    Thot_linpow = Ihot_long
                    Tcold_linpow = Icold_long
                    time_axis_raw_s = timeax
                    Tcold_freq = these_Tcold
                else:
                    raw_data = np.append(raw_data, Dlin, axis=1)
                    reduced_data = np.append(reduced_data, tmp_red_data, axis=1)
                    solar_obs = np.append(solar_obs, tmp_solar_obs, axis=1)
                    bkg_obs = np.append(bkg_obs, tmp_bkg_obs, axis=1)
                    sun_linpow = np.append(sun_linpow, sun, axis=1)
                    Thot_linpow = np.append(Ihot_long, Thot_linpow, axis=1)
                    Tcold_linpow = np.append(Icold_long, Tcold_linpow, axis=1)
                    time_axis_raw_s = np.append(time_axis_raw_s, new_time_axis)
                    Tcold_freq = np.append(Tcold_freq, these_Tcold)
                freq_axis_MHz      = freqs
                start_times        = np.append(start_times, time_start)
                time_axis_Thot_s   = np.append(time_axis_Thot_s, these_times_Thot)
                time_axis_Tcold_s  = np.append(time_axis_Tcold_s, these_times_Tcold)
                time_axis_sun_s    = np.append(time_axis_sun_s, these_times_sun)

                dt_time_axis = time_axis_raw_s[-1] + dT 


            ## Before storing the times, convert them from seconds to UTC
            time_axis_raw = pd.to_datetime(date+' '+start_times[0]) + pd.to_timedelta(time_axis_raw_s, unit='seconds') + pd.to_timedelta(dT, unit='seconds')
            time_axis_sun = pd.to_datetime(date+' '+start_times[0]) + pd.to_timedelta(time_axis_sun_s, unit='seconds') + pd.to_timedelta(dT, unit='seconds')
            time_axis_Thot = pd.to_datetime(date+' '+start_times[0]) + pd.to_timedelta(time_axis_Thot_s, unit='seconds') + pd.to_timedelta(dT, unit='seconds')
            time_axis_Tcold = pd.to_datetime(date+' '+start_times[0]) + pd.to_timedelta(time_axis_Tcold_s, unit='seconds') + pd.to_timedelta(dT, unit='seconds')

            ## Save the variables externally in a pkl file
            name_pkl = os.path.join(directory_output_day, filename_output_pkl) + '_' + date_stringformat + '.pkl'
            with open(name_pkl, 'wb') as pkl:
                pickle.dump(date, pkl)
                pickle.dump(raw_data, pkl)
                pickle.dump(reduced_data, pkl)
                pickle.dump(freq_axis_MHz, pkl)
                pickle.dump(time_axis_raw, pkl)
                pickle.dump(solar_obs, pkl)
                pickle.dump(bkg_obs, pkl)
                pickle.dump(Thot_linpow, pkl)
                pickle.dump(Tcold_linpow, pkl)
                pickle.dump(start_times, pkl)
                pickle.dump(time_axis_sun, pkl)
                pickle.dump(time_axis_Thot, pkl)
                pickle.dump(time_axis_Tcold, pkl)
                pickle.dump(Tcold_freq, pkl)
                pickle.dump(Thot, pkl)
        
        except Exception as e:
            
            print()
            print('Error message: ')
            print(e)
            print('!!!!!!!!!!!!!!!!!')
            print('An error occurred (see above) during the processing of the following day: '+this_date)
            print('Check if there was some telescope operation on that day (check the LogFile_telescope-operations.md file)')
            print('!!!!!!!!!!!!!!!!!')
            raise
            continue

    ## The following is just for debugging.
    # set which_plot = 0 to not doing anything

    which_plot = 0

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(dB[74,:], c='tab:blue')
    #ax.set_ylabel('dB')
    #axr = ax.twinx()
    #axr.plot(data[74,:], c='tab:orange', linestyle='dashed')
    #axr.set_ylabel('raw dates', c='tab:orange')
    #plt.show()
    
    if which_plot == 1:
            
        #idx_f1 = np.abs(freq_axis_MHz - 300).argmin()
        #idx_f2 = np.abs(freq_axis_MHz - 600).argmin()
        idx_f1 = np.abs(freq_axis_MHz - 11060).argmin()
        idx_f2 = np.abs(freq_axis_MHz - 10640).argmin()
        df = 5
        dfmin1 = idx_f1 - df
        dfmax1 = idx_f1 + df
        dfmin2 = idx_f2 - df
        dfmax2 = idx_f2 + df
        lc_sun = np.mean(sun_linpow[dfmin1:dfmax1,:],axis=0)
        lc_sun2 = np.mean(sun_linpow[dfmin2:dfmax2,:],axis=0)
        lc_Thot = np.mean(Thot_linpow[dfmin1:dfmax1,:],axis=0)
        lc_Thot2 = np.mean(Thot_linpow[dfmin2:dfmax2,:],axis=0)
        lc_Tcold = np.mean(Tcold_linpow[dfmin1:dfmax1,:],axis=0)
        lc_Tcold2 = np.mean(Tcold_linpow[dfmin2:dfmax2,:],axis=0)

        Yfact1 = lc_Thot / lc_Tcold
        Yfact1 = 10.0*np.log10(Yfact1)
        Yfact2 = lc_Thot2 / lc_Tcold2
        Yfact2 = 10.0*np.log10(Yfact2)

        fig = plt.figure(figsize=(7,4))

        ax1 = fig.add_subplot(211)
        #extent = time_axis_raw_s[0], time_axis_raw_s[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
        #im = ax1.imshow(raw_data, aspect='auto', extent=extent)
        #cbar = fig.colorbar(im, ax=ax1)
        #cbar.set_label('Linear power')
        ax1.set_title('freq1')
        ax1.plot(time_axis_sun, lc_sun, '.', c='orange', alpha=0.8, label='Sun')
        ax1.plot(time_axis_Tcold, lc_Tcold, '.', c='blue', alpha=0.8, label='Tcold')
        ax1.plot(time_axis_Thot, lc_Thot, '.', c='red', alpha=0.8, label='Thot')
        #ax1.set_xlabel('Time')
        ax1.set_ylabel('linear power')
        ax1.legend()
        #ax1r = ax1.twinx()
        #ax1r.plot(time_axis_Tcold_s, Yfact1, '.', c='gray')
        #ax1r.set_ylabel('Y-factor [-]')
        #ax1r.set_yscale('log')

        ax2 = fig.add_subplot(212)
        #ax2.plot(time_axis_sun_s, lc_sun, '.', c='tab:orange', label='11.06 GHz')
        #ax2.plot(time_axis_sun_s, lc_sun2, '.', c='tab:green', label='10.6 GHz')
        ax2.set_title('freq2')
        ax2.plot(time_axis_sun, lc_sun2, '.', c='orange', alpha=0.8, label='Sun')
        ax2.plot(time_axis_Tcold, lc_Tcold2, '.', c='blue', alpha=0.8, label='Tcold')
        ax2.plot(time_axis_Thot, lc_Thot2, '.', c='red', alpha=0.8, label='Thot')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Linear power')
        ax2.legend()
        #ax2r = ax2.twinx()
        #ax2r.plot(time_axis_Tcold_s, Yfact2, '.', c='gray')
        #ax2r.set_ylabel('Y-factor [-]')
        #ax2r.set_yscale('log')

        plt.show()
    elif which_plot == 2:

        fig = plt.figure(figsize=(7,4))

        ax1 = fig.add_subplot(211)
        #extent = timeax[0], timeax[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
        #im = ax1.imshow(tmp_solar_obs, aspect='auto', vmin=100, vmax=600, extent=extent)
        extent = time_axis_raw[0], time_axis_raw[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
        im = ax1.imshow(solar_obs, aspect='auto', vmin=100, vmax=600, extent=extent)
        cbar = fig.colorbar(im, ax=ax1)
        cbar.set_label('Solar flux [sfu]')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Freq [GHz]')
        
        ax2 = fig.add_subplot(212)
        #extent = timeax[0], timeax[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
        #im = ax1.imshow(tmp_solar_obs, aspect='auto', vmin=100, vmax=600, extent=extent)
        extent = time_axis_raw[0], time_axis_raw[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
        im2 = ax2.imshow(raw_data, aspect='auto', extent=extent, vmax=200)
        cbar = fig.colorbar(im2, ax=ax2)
        cbar.set_label('Linear power')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Freq [GHz]')

        plt.show()

if __name__ == '__main__':
    main()