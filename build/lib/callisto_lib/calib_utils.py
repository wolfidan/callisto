import numpy as np
import pandas as pd

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

