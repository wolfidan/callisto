'''
Functions to get the estimated Tcold. As the determination of the cloudiness of the sky
is made by humans, this quantity cannot be downloaded live. The table in get_clouds needs
to be updated regularly. To do so, follow these steps:

    1. From CLIMAP, open the data 'Gesamtbewölkung' (nto000s0) from Locarno-Monti as table.
    2. Select all columns and rows and copy the data.
      !! It is important that only time and Gesamtbewölkung are copied !!
    3. Paste the data in the following file:
        'C:\\xrt\\calibration\\Tcold\\measured_cloudiness.txt'
    (the old content can be deleted)
    4. run update_csv_clouds()

Author: Andrea F. Battaglia (2025-01-20)
'''

import pandas as pd
import numpy as np


### The following code has been commented out because it is not used anymore
'''
def update_csv_clouds():

    ### Path to the CSV containing the octas of clouds
    path_csv_clouds = 'C:\\xrt\\calibration\\Tcold\\clouds_octas.csv'

    ### Data obtained from CLIMAP
    in_path_climap = 'C:\\xrt\\calibration\\Tcold\\measured_cloudiness.txt'

    ### Read the data
    df_climap = pd.read_csv(in_path_climap, sep='\t', header=None)
    times_climap = df_climap[0].values
    val_climap = df_climap[1].values
    times_climap = times_climap[val_climap >= 0]
    times_climap = [pd.to_datetime(str(t), format='%d.%m.%Y %H:%M:%S') for t in times_climap]
    val_climap = val_climap[val_climap >= 0]

    ### Create the new CSV file if it does not exist, otherwise append
    df_climap_new = pd.DataFrame({'times': times_climap, 'octas': val_climap})
    df_climap_new = df_climap_new.sort_values(by='times')
    df_climap_new = df_climap_new.reset_index(drop=True)
    try: df_climap_new.to_csv(path_csv_clouds, mode='x', index=False)
    except: df_climap_new.to_csv(path_csv_clouds, mode='a', header=False, index=False)

    ### Remove any duplicates from the file
    df_climap_new = pd.read_csv(path_csv_clouds)
    df_climap_new = df_climap_new.drop_duplicates()
    df_climap_new.to_csv(path_csv_clouds, index=False)


def get_clouds(time, print_clouds=False):

    ### Path to the CSV containing the octas of clouds
    path_csv_clouds = 'C:\\xrt\\calibration\\Tcold\\clouds_octas.csv'

    ### Read the CSV file
    df_clouds = pd.read_csv(path_csv_clouds)
    times_clouds = df_clouds['times'].values
    times_clouds = [pd.to_datetime(t) for t in times_clouds]
    octas = df_clouds['octas'].values

    ### Get the closest time
    diff = [abs(t - time) for t in times_clouds]
    idx = np.argmin(diff)
    time_closest = times_clouds[idx]

    if print_clouds: print(f'The closest time is {time_closest} with {octas[idx]} octas of clouds')

    return octas[idx], time_closest
'''


def get_Tcold(time, myFreq, myElev, precip_mm, modeln=1):
    '''
    Get the estimated Tcold, by interpolating the values from the CSV file

    Parameters
    ----------
    time : datetime or str
        Time of the observation.

    myFreq : array_like
        Frequency of the observation in GHz.

    myElev : float
        Elevation of the observation in degrees.
        
    precip_mm : float
        Precipitation [mm] near the time of interest.

    modeln : int, optional
        Model number (for the CSV). The default is 1.
            - model 1: AH = 15 g/m3, Columnar Liquid = 0.5 kg/m2
            - model 2: AH = 20 g/m3, Columnar Liquid = 1.2 kg/m2
    '''

    time = pd.to_datetime(time)
    
    ### CSV file containing the tabulated values
    path_csv = 'C:\\xrt\\calibration\\Tcold\\Tcold_look_up_table.csv'
    
    ### First of all, we upèdate the table of clouds, in the case new data came in
    #update_csv_clouds()
    
    ### Read Tcold from the CSV file
    df_Tcold = pd.read_csv(path_csv, header=15)

    ### Select the desired model
    if modeln == 1: df_Tcold_model = df_Tcold[:28] 
    elif modeln == 2: df_Tcold_model = df_Tcold[28:]

    ### Apply the criterion to determine clear or cloudy sky
    #****************************************************************************************************
    ## This is a temporary way we used, based on octas determination by observers. As this criterion
    ## can be a bit subjective, we abandoned it. Now, we only say that it is not clear sky based on the
    ## measured precipitation (see below).
    #eps_clouds = 4  # threshold clear/cloudy sky in octas
    #these_coulds, time_clouds = get_clouds(time)
    #if time.date() != time_clouds.date():
    #    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    #    print('WARNING! The dates of the CALLISTO and clouds data are different.')
    #    print('Probably, updated CLIMAP data needs to be downloaded.')
    #    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    #if these_coulds <= eps_clouds: df_Tcold_model = df_Tcold_model[df_Tcold_model['cond'] == 'clear']
    #else: df_Tcold_model = df_Tcold_model[df_Tcold_model['cond'] == 'cloudy']
    #****************************************************************************************************
    ## This is a temporary way we used, based on octas determination by observers. As this criterion
    ## can be a bit subjective, we abandoned it. Now, we only say that it is not clear sky based on the
    ## measured precipitation (see below).
    eps_precip = 0.  # threshold clear/cloudy sky in mm of precipitation
    if precip_mm <= eps_precip: df_Tcold_model = df_Tcold_model[df_Tcold_model['cond'] == 'clear']
    else: df_Tcold_model = df_Tcold_model[df_Tcold_model['cond'] == 'cloudy']    
    
    ### Do linear interpolation
    Tcold_10GHz = df_Tcold_model[df_Tcold_model['freq'] == 10]['Tcold'].values[::-1]
    elev_10GHz = df_Tcold_model[df_Tcold_model['freq'] == 10]['elev'].values[::-1]
    Tcold_12GHz = df_Tcold_model[df_Tcold_model['freq'] == 12]['Tcold'].values[::-1]
    elev_12GHz = df_Tcold_model[df_Tcold_model['freq'] == 12]['elev'].values[::-1]  

    elev = np.linspace(0, 90, 100)
    Tcold_10GHz_interp = np.interp(elev, elev_10GHz, Tcold_10GHz)
    Tcold_12GHz_interp = np.interp(elev, elev_12GHz, Tcold_12GHz)

    ### Get the index of the closest elevation to myElev for both
    idx_10GHz = np.argmin(np.abs(elev - myElev))
    idx_12GHz = np.argmin(np.abs(elev - myElev))
    Tcold_10GHz_interp[idx_10GHz], Tcold_12GHz_interp[idx_12GHz]

    ### linearly interpolate to myFreq with np.interp
    Tcold_myFreq_interp = [np.interp(myFreq[ff], [10, 12], [Tcold_10GHz_interp[idx_10GHz], Tcold_12GHz_interp[idx_12GHz]]) for ff in range(len(myFreq))]

    return Tcold_myFreq_interp

