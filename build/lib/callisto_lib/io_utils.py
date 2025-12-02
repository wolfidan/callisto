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


def read_motor_log(day, path_log_motor = None):
    '''
    This function is used to read the motor DISEQ log for a given day
    '''
    if not path_log_motor:
        path_log_motor = os.path.join(constants.PATH_LOG, "DISEQ")
        
    date_motor = day[0:4]+'-'+day[4:6]+'-'+day[6:]
    filename_format_log = f'DISEQ-{date_motor}-Sun.txt'
    path_log_motor = os.path.join(path_log_motor, filename_format_log)
    df_motor = pd.read_csv(path_log_motor, sep=',', header=3)
    df_motor.rename(columns=lambda x: x.strip(), inplace=True)
    return df_motor


#**************************************************************************

def get_meteoswiss_data(date, path_meteoswiss_data = None):
    """
    Download the OGDS-SMN 't_now' CSV file for the given station and 
    return it as a pandas DataFrame.

    Parameters
    ----------
    station : str
        Station name (case-insensitive)

    Returns
    -------
    pd.DataFrame
    """
    station_l = station.lower()
    url = f"https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn/{station_l}/ogd-smn_{station_l}_t_now.csv"

    resp = requests.get(url)
    resp.raise_for_status()  # raise error if download fails

    # Read CSV from memory
    csv_buffer = StringIO(resp.text)
    df = pd.read_csv(csv_buffer)

    return df
        
#**************************************************************************

def get_fit_files(day, path_fit_folder = None):
    '''
    This function is used to get a list of all fit files and their times
    for a given day
    '''
    if not path_fit_folder:
        path_fit_folder = os.path.join(constants.PATH_DATA, "raw", "FITfiles")
        
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

def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]