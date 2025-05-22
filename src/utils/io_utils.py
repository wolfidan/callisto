import datetime
from astropy.io import fits
import numpy as np
import os
import pandas as pd
import glob

# Local imports
import constants

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
    
    # add datetime for convenience
    date = datetime.datetime.strptime(dict_fitfile["date"] 
            + "_" + dict_fitfile["T0"], "%Y/%m/%d_%H:%M:%S.%f")
    dict_fitfile["datetimes"] = np.array([(date + 
                datetime.timedelta(seconds = dt)).replace(tzinfo=datetime.UTC)
                                 for dt in dict_fitfile["timeax"]])
    return dict_fitfile


def log_execution(fname, message, plot_separation=True):
    """
    Description:
        This is to store in an external file a log, whether the script run
        successfully or an error occurred.
    """
    with open(fname, 'a') as log_file:
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_file.write(f'\n{timestamp} - {message}')
        if plot_separation == True: 
            log_file.write(f'\n************************************************\n')

def read_motor_log(day, path_log_motor = constants.PATH_LOG_MOTOR):
    '''
    This function is used to read the motor DISEQ log for a given day
    '''
    date_motor = day[0:4]+'-'+day[4:6]+'-'+day[6:]
    filename_format_log = f'DISEQ-{date_motor}-Sun.txt'
    path_log_motor = os.path.join(path_log_motor, filename_format_log)
    df_motor = pd.read_csv(path_log_motor, sep=',', header=3)
    df_motor.rename(columns=lambda x: x.strip(), inplace=True)
    return df_motor


#**************************************************************************

def get_meteoswiss_data(date):
    ### Find the meteoswiss CSV file and open it
    filename_csv = f'meteoswiss_{date.strftime("%Y-%m-%d")}.csv'
    filepath_csv = os.path.join(constants.PATH_METEOSWISS_DATA, filename_csv)
    df_meteoswiss = pd.read_csv(filepath_csv, sep=',')
    df_meteoswiss['time'] = pd.to_datetime(df_meteoswiss['time'])
    
    return df_meteoswiss
        
#**************************************************************************

def get_fit_files(day, path_fit_folder = constants.PATH_FIT_FOLDER):
    '''
    This function is used to get a list of all fit files and their times
    for a given day
    '''
    filename_format = 'meteoswiss_{:s}_01.fit' #MyFile = 'SWISS-METEO_20231127_*_01.fit'
    filepath_format = os.path.join(path_fit_folder, filename_format.format(day + "_*"))
    list_fit = sorted(glob.glob(filepath_format))
    list_fit.sort()
    # Get time of all fit files
    list_fit_times = [datetime.datetime.strptime(os.path.basename(fitname), 
                    filename_format.format("%Y%m%d_%H%M%S")).replace(tzinfo=datetime.UTC)
                    for fitname in list_fit]
    list_fit = np.array(list_fit)
    list_fit_times = np.array(list_fit_times)
    return list_fit, list_fit_times
