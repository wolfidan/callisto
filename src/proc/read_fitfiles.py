# -*- coding: utf-8 -*-

description = """
Description:

    Read FIT-files from local folder.

    @author: Philipp Schmid

    2024/05/10: make draft based on Christian Monstein's script "PlottSunLCfromFIT-loop.py"

    #### compute Aeff at each frequency
    #### define frequencies over which to integrate
    #### ...
    
    History:
        - 2024/12/11 [Andrea F. Battaglia]: add some comments for clarity
        - 2024/12/12 [Andrea F. Battaglia]: now the script does not process all raw data. It
                                            checks which data has no counterparts in the
                                            output directory of the processed observations.
                                            However, if for whathever reason all files have
                                            to be reprocessed, run "-day all"
        - 2024/12/13 [Andrea F. Battaglia]: in addition to the standard csv file, now the 
                                            software now stores externally a pkl file
                                            containing some useful variables that allow to
                                            look at the raw data at high-resolution. This is
                                            helpful to investigate systematic behaviours.
                                            Also, the Thot, Tcold and solar measurements are
                                            stored in separate variables
        - 2024/12/17 [Andrea F. Battaglia]: Bug fix (it had issues in generating the pkl files
                                            when the argument '-raw all' was run)

"""

import numpy as np
import scipy.constants as CON
import argparse
import glob
import sys
import os
import matplotlib.pyplot as plt
import pickle

from tools.return_fits2process import *
from tools.utils_fitfile import readfit
from tools.utils import checkdirs

def compute_Aeff(gain_dB, freq):

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
    
    parser.add_argument(
        '-verbose', action='store_true', required=False, default=False,
        help='Print more information like results (default: %(default))')
    
    parser.add_argument(
        '-day', type=str, metavar=str.__name__, 
        action='store', required=False, default=None,
        help='Choose day to analyze (default: analyze all days)')
    
    args = parser.parse_args()
    verbose = args.verbose
    day = args.day

    # --------------------------
    # Definitions (to be revised and adjusted!)
    # --------------------------

    Tcold = 12 # from literature
    Thot = 273.15 + 11.0 # from outdoor sensor
    Terr = 0 # timing error FIT-files in seconds, should be 0 in normal cases
    
    gain_dB = 36 # according to data sheet of TRIAX TDS 65 A
    sfu = 1e22

    freq_MHz = 11075 # in case of 0.25 sec FIT-file
    dfreq_MHz = 50 # MHz
    # [0081]=11056.562,0
    # [0082]=11064.188,0
    # [0083]=11071.812,0
    # [0084]=11079.438,0
    # [0085]=11087.062,0
    # [0086]=11094.688,0

    ### that should be done for each frequency, no?!!!?!!
    Aeff = compute_Aeff(gain_dB, freq_MHz*1e6) 
    gain_linear = 10.0**(gain_dB/10.) 
    wavelength = CON.c/(freq_MHz*1e6)
    Aeff = wavelength**2 * gain_linear / (4*np.pi)
    
    # --------------------------
    # Define files to be read
    # Create output directory 
    # --------------------------

    filename_output = 'solarflux.csv'
    filename_output_pkl = 'time-profiles_sun_Thot_Tcold'
    directory_fitfiles = 'C:\\xrt\\output\\data\\raw\\FITfiles'
    directory_output = 'C:\\xrt\\output\\data\\proc'
    checkdirs([directory_fitfiles, directory_output])

    ## If the user run the script with the argument "-day" activated, then enter this
    if day:
        if day == 'all':
            ## If "-day all" is run, then reprocess all files
            filename_format = 'meteoswiss_*_*_01.fit'
        else:
            ## else, only the selected day
            filename_format = 'meteoswiss_%s_*_01.fit' % day #MyFile = 'SWISS-METEO_20231127_*_01.fit'
        filepath_format = os.path.join(directory_fitfiles, filename_format)
        filepaths = sorted(glob.glob(filepath_format))
    else:
        ## Detect which files have not been processed yet
        files2process = return_fits2process(directory_fitfiles, directory_output)
        filepaths = [os.path.join(directory_fitfiles, f) for f in files2process]
        filepaths.sort()

    
    if len(filepaths) < 1:
        print('All raw data has been processed (or no FIT-file has been found). Return without processing')
        return
        #print('No FIT-files found: %s' % filepath_format )
        #sys.exit()
    # else:
    #     if verbose:
    #         print('FIT-files found: ',(len(filepaths)),filepaths)

    ## Create the folder where to store the processed data. In addition
    # create the csv file where some information will be stored.
    # Here we also initialize the csv file by defining its header, containing
    # some calibration information
    dates = list(set([os.path.basename(path)[11:19] for path in filepaths]))
    for date in dates:
        directory_output_day = os.path.join(directory_output, date)
        if not os.path.isdir(directory_output_day):
            print('Create directory: %s' % directory_output_day)
            os.makedirs(directory_output_day)
        print('Create file: %s' % os.path.join(directory_output_day, filename_output))

        # calibration information
        with open(os.path.join(directory_output_day, filename_output), 'w') as f:
            f.write('# Tcold_K: %s\n' % Tcold)
            f.write('# Thot_K: %s\n' % Thot)
            f.write('# freq_MHz: %s\n' % freq_MHz)
            f.write('# dfreq_MHz: %s\n' % dfreq_MHz)
        
        # results
        header = ['date', 'time_start',  
                  'Icold_adu', 'Ihot_adu', 'Ycal_dB', 
                  'Isun_mean_adu', 'Isun_max_adu', 'Isun_mean_above_max_adu', 'nrvals_for_Isun_mean_above_max',
                  'Ssun_sfu', 'Sback_sfu', 'Ssun_minus_Sback_sfu']
        header_txt = ', '.join(header) + '\n'
        with open(os.path.join(directory_output_day, filename_output), 'a') as f:
            f.write(header_txt)

    # --------------------------
    # Process files
    # --------------------------

    ## We initialize some variables that then will be stored in a separate pkl file.
    # These variables are useful to diagnose the observations and assess the calibration
    # status during a particular day.
    raw_data           = []
    freq_axis_MHz      = []
    time_axis_raw_s    = []
    start_times        = []
    lc_sun_linpow      = []
    lc_Thot_linpow     = []
    lc_Tcold_linpow    = []
    time_axis_sun_s    = []
    time_axis_Thot_s   = []
    time_axis_Tcold_s  = []

    dict_fitfile_tmp = readfit(filepaths[0])
    date = dict_fitfile_tmp['date'] # e.g. 2024/05/10
    old_date_stringformat = date[0:4] + date[5:7] + date[8:] # 20240510
    enter_now = True
    ## Loop over all files (one file is generated every 15 minutes)
    #for filepath in filepaths:
    for ii in range(len(filepaths)):
        filepath = filepaths[ii]

        ## Here we define the freqeuncy band over which to integrate
        # It would be good to include the processing of different freqeuncy bands
        # (hence, to loop also over frequencies)
        # indices 81-86 correspond to the output frequency of 11.075 GHz    
        freqidx_min = 81
        freqidx_max = 86

        ## Read the fit files
        dict_fitfile = readfit(filepath)
        data = dict_fitfile['data']         # array of values in raw-date
        freqs = dict_fitfile['freqs']       # (original) requency axis in GHz
        timeax = dict_fitfile['timeax']     # time axis in s
        time_start = dict_fitfile['T0']
        dT = dict_fitfile['dT']
        date = dict_fitfile['date'] # e.g. 2024/05/10
        date_stringformat = date[0:4] + date[5:7] + date[8:] # 20240510
        directory_output_day = os.path.join(directory_output, date_stringformat)

        if date_stringformat != old_date_stringformat:
            raw_data           = []
            freq_axis_MHz      = []
            time_axis_raw_s    = []
            start_times        = []
            lc_sun_linpow      = []
            lc_Thot_linpow     = []
            lc_Tcold_linpow    = []
            time_axis_sun_s    = []
            time_axis_Thot_s   = []
            time_axis_Tcold_s  = []

        ## This is to make the time axis increasing, to easily plot the time profiles
        if date_stringformat != old_date_stringformat or enter_now:
        #if filepath == filepaths[0]:
            new_time_axis = timeax
        else: 
            new_time_axis = timeax+time_axis_raw_s[-1]+dT

        ## Convert the raw-dates to dB and then to linear power
        dB = data/255*2500.0/25.4 # raw-date -> dB
        Dlin = 10.0**(dB/10) # dB -> linear power
        #T = np.arange(0,len(Dlin[0]),1) * dT # sampling rate usually 2 Hz
        
        ## Average the linear power array over a freqency band
        # LC = np.mean(Dlin[140:180,:], axis=0) # take only frequencies without SAT-TV / 0.5 sec FIT-files with 200 channels
        # LC = np.mean(Dlin[40:80,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files with 100 channels
        LC = np.mean(Dlin[freqidx_min:freqidx_max,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files with 100 channels
        
        ## cold calibration
        T1 = int((10-Terr)/dT) # start Icold, sec -> pixel
        T2 = int((60-Terr)/dT)
        Icold = np.mean(LC[T1:T2])
        Icold_long = LC[T1:T2]
        these_times_Tcold = new_time_axis[T1:T2]

        ## hot calibration
        T3 = int(( 80-Terr)/dT) # start Ihot, sec -> pixel
        T4 = int((120-Terr)/dT)
        Ihot  = np.mean(LC[T3:T4])
        Ihot_long = LC[T3:T4]
        these_times_Thot = new_time_axis[T3:T4]
        
        ## Y value
        Ycal  = Ihot/Icold
        YcaldB = 10.0*np.log10(Ycal)

        # # calibration snr
        # Icold_std = np.std(LC[T1:T2])
        # snr = (Ihot-Icold)/Icold_std
        # snr_dB = 10*np.log10(snr)
        
        ## Extract the sun measurements
        T5 = int((170-Terr)/dT) # start Isunscan, sec -> pixel
        T6 = int((900-Terr)/dT) # 900
        sun = LC[T5:T6]
        these_times_sun = new_time_axis[T5:T6]
               
        MyFluxLimit = 0.97 # take only values above 97% of peak flux
        if len(sun) > 0:
            Isun_mean = np.mean(sun)
            Isun_max = np.max(sun)
            Isun_mean_above_max = np.mean(sun[sun>MyFluxLimit*Isun_max])
            nrvals_for_Isun_mean_above_max = np.sum(sun>MyFluxLimit*Isun_max)
        else:
            Isun_mean = np.nan
            Isun_max = np.nan
            Isun_mean_above_max = np.nan
            nrvals_for_Isun_mean_above_max = 0
      
        Tsun = (Thot - Tcold) / (Ihot - Icold) * (Isun_mean_above_max - Icold) + Tcold # Antenna temperature in Kelvin
        Ssun = 2*CON.Boltzmann*Tsun/Aeff*sfu
        Sback = (2*CON.Boltzmann*Tcold)/(wavelength**2*gain_linear)*4*np.pi * sfu
        Ssun_minus_Sback = Ssun - Sback
        
        #-------------------------------------------------------------------------------
        
        ## Store the variables in loop
        if date_stringformat != old_date_stringformat or enter_now:
        #if filepath == filepaths[0]:
            raw_data = data
            time_axis_raw_s = timeax
            enter_now = False
        else:
            raw_data = np.append(raw_data, data, axis=1)
            time_axis_raw_s = np.append(time_axis_raw_s, new_time_axis)
        freq_axis_MHz = freqs
        start_times = np.append(start_times, time_start)
        lc_sun_linpow  = np.append(lc_sun_linpow, sun)
        lc_Thot_linpow     = np.append(lc_Thot_linpow, Ihot_long)
        lc_Tcold_linpow    = np.append(lc_Tcold_linpow, Icold_long)
        time_axis_Thot_s   = np.append(time_axis_Thot_s, these_times_Thot)
        time_axis_Tcold_s  = np.append(time_axis_Tcold_s, these_times_Tcold)
        time_axis_sun_s    = np.append(time_axis_sun_s, these_times_sun)

        if verbose:
            print('Icold {:8.1f} digit, Ihot {:8.1f} digit, YcaldB {:5.2f} dB '.format(Icold,Ihot,YcaldB))
            print('Solar radio flux:            {:6.1f} SFU'.format(Ssun))
            print('Background flux:           {:6.1f} SFU'.format(Sback))
            print('Final flux above background: {:6.1f} SFU'.format(Ssun_minus_Sback))
        
        newline = [date, time_start, Icold, Ihot, YcaldB, 
                   Isun_mean, Isun_max, Isun_mean_above_max, nrvals_for_Isun_mean_above_max,
                   Ssun, Sback, Ssun_minus_Sback]
        newline_txt = ', '.join([str(el) for el in newline]) + '\n'
        with open(os.path.join(directory_output_day, filename_output), 'a') as f:
            f.write(newline_txt)

        if filepath != filepaths[-1]: 
            next_file = readfit(filepaths[ii+1])
            date_next_file = next_file['date']
            next_date_stringformat = date_next_file[0:4] + date_next_file[5:7] + date_next_file[8:] # 20240510
            
            ## If the next file is on a different date, store the pkl file of today
            if next_date_stringformat != date_stringformat:
                
                ## Save the variables externally in a pkl file
                name_pkl = os.path.join(directory_output_day, filename_output_pkl) + '_' + date_stringformat + '.pkl'
                with open(name_pkl, 'wb') as pkl:
                    pickle.dump(date, pkl)
                    pickle.dump(raw_data, pkl)
                    pickle.dump(freq_axis_MHz, pkl)
                    pickle.dump(time_axis_raw_s, pkl)
                    pickle.dump(start_times, pkl)
                    pickle.dump(lc_sun_linpow, pkl)
                    pickle.dump(lc_Thot_linpow, pkl)
                    pickle.dump(lc_Tcold_linpow, pkl)
                    pickle.dump(time_axis_sun_s, pkl)
                    pickle.dump(time_axis_Thot_s, pkl)
                    pickle.dump(time_axis_Tcold_s, pkl)
        else:
            ## Save the variables externally in a pkl file
            name_pkl = os.path.join(directory_output_day, filename_output_pkl) + '_' + date_stringformat + '.pkl'
            with open(name_pkl, 'wb') as pkl:
                pickle.dump(date, pkl)
                pickle.dump(raw_data, pkl)
                pickle.dump(freq_axis_MHz, pkl)
                pickle.dump(time_axis_raw_s, pkl)
                pickle.dump(start_times, pkl)
                pickle.dump(lc_sun_linpow, pkl)
                pickle.dump(lc_Thot_linpow, pkl)
                pickle.dump(lc_Tcold_linpow, pkl)
                pickle.dump(time_axis_sun_s, pkl)
                pickle.dump(time_axis_Thot_s, pkl)
                pickle.dump(time_axis_Tcold_s, pkl)

        old_date_stringformat = date_stringformat

    #plt.figure()
    #plt.plot(time_axis_sun_s, lc_sun_linpow,'.',markersize=1)
    #plt.plot(time_axis_Tcold_s, lc_Tcold_linpow,'.',markersize=1)
    #plt.plot(time_axis_Thot_s, lc_Thot_linpow,'.',markersize=1)
    #plt.show()

if __name__ == '__main__':
    main()
    