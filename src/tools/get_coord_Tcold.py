'''
Get the coordinates of Tcold for a given measurement

Author: Andrea F. Battaglia (2025-01-20)
'''

import pandas as pd
import os
import numpy as np


def get_DISEQ_file(datetime):
    '''
    Get the DISEQ file for a given time

    Parameters
    ----------
    datetime : str or pandas.Timestamp
        Time of the measurement

    Returns
    -------
    str
        Path to the DISEQ file
    '''

    ### Folder containing the files with the rotor coordinates
    folder_rotor_coord = 'C:\\xrt\\src\\PythonScripts\\TrackingSun\\'

    ### Find all DISEQ files
    files = os.listdir(folder_rotor_coord)
    files = [f for f in files if 'DISEQ-' in f and 'Sun.txt' in f]
    files.sort()

    ### Get the dates of the files    
    date_files = [f.split('-')[1]+'-'+f.split('-')[2]+'-'+f.split('-')[3] for f in files]

    ### Find the file
    idx = date_files.index(datetime.date().strftime('%Y-%m-%d'))
    filename = 'DISEQ-'+date_files[idx]+'-Sun.txt'

    return folder_rotor_coord+filename


def get_AziElev_Tcold(datetime):
    '''
    Get the azimuth and elevation of Tcold for a given measurement

    Parameters
    ----------
    datetime : str or pandas.Timestamp
        Time of the measurement

    Returns
    -------
    float, float
        Azimuth and elevation of Tcold
    '''

    ### Get the DISEQ file
    datetime = pd.to_datetime(datetime)
    file_diseq = get_DISEQ_file(datetime)

    ### Read the file
    data_motor = np.genfromtxt(file_diseq, dtype=str, delimiter=',', 
                    skip_header=5, filling_values='', usecols=(0, 1, 2, 3),
                    encoding='latin-1')

    ### Get the data arrays
    time = [float(data_motor[i,0]) for i in range(len(data_motor[:,0]))]
    azi  = [float(data_motor[i,1]) for i in range(len(data_motor[:,1]))]
    ele  = [float(data_motor[i,2]) for i in range(len(data_motor[:,2]))]
    comm = [data_motor[i,3] for i in range(len(data_motor[:,3]))]

    ### Select only the indices where comm contains 'Cold sky reference' in the whole string
    idx = [i for i in range(len(comm)) if 'Cold sky reference' in comm[i]]
    time = [time[i] for i in idx]
    azi = [azi[i] for i in idx]
    ele = [ele[i] for i in idx]

    ### COnvert the hour times to pandas.Timestamp
    hhmm = []
    for t in time:
        h = int(t)
        m = int((t-h)*60)
        hhmm.append('{:02d}:{:02d}'.format(h,m))
    dates_motor = [pd.to_datetime(str(datetime.date()) + ' ' + t) for t in hhmm]

    ### Get the closest time
    diff = [abs(dates_motor[i] - datetime) for i in range(len(dates_motor))]
    idx = diff.index(min(diff))

    return azi[idx], ele[idx]