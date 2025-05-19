
import numpy as np
import datetime
from astropy.io import fits
import pandas as pd
import os
import glob

### Path to the folder where to get the log of the motors
PATH_LOG_MOTOR = 'C:\\xrt\\src\\PythonScripts\\TrackingSun'
### Path to the folder where to get the FIT files
PATH_FIT_FOLDER = 'C:\\xrt\\output\\data\\raw\\FITfiles'

#**************************************************************************

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

#**************************************************************************

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
    dict_fitfile["datetimes"] = np.array([date + datetime.timedelta(seconds = dt) 
                                 for dt in dict_fitfile["timeax"]])
    return dict_fitfile

#**************************************************************************

def convert_decimal_hour_timedelta(time_df):
    '''
    This function converts the times of the log file of the motors to timedelta
    '''
    hhmmss = []
    for t in time_df:
        
        # Extract hours
        hh = int(t)
        # Extract minutes
        mm = int((t - hh) * 60)
        # Extract seconds
        ss = int(((t - hh) * 60 - mm) * 60)
        # Extract milliseconds
        mil = int(((t - hh) * 60 - mm - ss/60) * 60000)
        # Format time 
        hhmmss.append(datetime.timedelta(hours = hh, minutes = mm, seconds = ss,
                                         milliseconds = mil))
    
    return np.array(hhmmss)

#**************************************************************************

def gaussian_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D Gaussian to the solar scan image.
    This function is now (2025-03-31) obsolete. We now use the 2D paraboloid.
    '''
    x, y = coords
    return amplitude * np.exp(-(((x-x0)**2)/(2*sigma_x**2) + ((y-y0)**2)/(2*sigma_y**2)))

#**************************************************************************

def paraboloid_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D paraboloid to the solar scan image.
    '''
    x, y = coords
    return amplitude - 4*np.log(2) * ((x-x0)**2/sigma_x**2 + (y-y0)**2/sigma_y**2)

#**************************************************************************

def read_motor_log(day, path_log_motor = PATH_LOG_MOTOR):
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

def get_fit_files(day, path_fit_folder = PATH_FIT_FOLDER):
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
                    filename_format.format("%Y%m%d_%H%M%S"))
                    for fitname in list_fit]
    return list_fit, list_fit_times

#**************************************************************************

def average_around_time(times, data, ref_time, avg_period = 1.5):
    dic = {"timestamp": times, "value": data}
    df = pd.DataFrame(dic)
    
    # Ensure timestamp is a datetime type and sorted
    df['timestamp'] = pd.to_datetime(df['timestamp'])
    df = df.sort_values('timestamp')

    # Set timestamp as index for efficient slicing
    df = df.set_index('timestamp')

    # Duration of the averaging window
    delta = pd.Timedelta(seconds=avg_period)

    # Function to compute mean in the time window around a given time
    def average_around_date(date):
        window = df.loc[date - delta : date + delta]
        return window['value'].mean()

    # Apply the function to each date
    averages = [average_around_date(i) for i in ref_time]

    return np.array(averages)

