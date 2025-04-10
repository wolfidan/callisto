'''
This script, which is run every 15 minutes, allows to get in real-time (from the Locarno-Monti stations):
   - temperature [degC] 2 meters above ground, 
   - the precipitation [mm], and 
   - the global radiation (e.g. solar irradiance) [W/m2]. 
The json files are updated every 10 minutes (5 min delay between measure and update of the file).
Afterwards, these quantities are added to a csv file (daily). This allows then to
process every FIT file taking into account proper weather conditions.

In data['features'][XX], the number XX stands for the station. The Locarno Monti stations are
XX=78 for the temperature and solar radiation stations and XX=77 for the precipitation station.

Created on 2025-01-21 by Andrea F. Battaglia
'''

import sys
sys.stderr = open('error_log_meteoswiss_data.txt', 'w')

import urllib.request
import json
import pandas as pd
import os
import datetime

from astropy.io import fits


def create_csv_file(csv_filename):
    with open(csv_filename, 'w') as f:
        f.write('time,temp_degC,precip_mm,irrad_W_m2\n')
    return csv_filename


def get_T_Locarno_Monti(print_T=False):
    json_link_T = 'https://data.geo.admin.ch/ch.meteoschweiz.messwerte-lufttemperatur-10min/ch.meteoschweiz.messwerte-lufttemperatur-10min_en.json'
    station_index = 78 # Locarno Monti
    
    # Get the JSON file
    with urllib.request.urlopen(json_link_T) as url:
        data = json.loads(url.read().decode())

    # Print, if requested
    if print_T:
        print('Station: '+data['features'][station_index]['properties']['station_name'])
        print('Time of the measurement: '+data['features'][station_index]['properties']['reference_ts'])
        print('Temperature at 2 meters above ground: '+str(data['features'][station_index]['properties']['value'])+' °C')

    # Return the temperature
    return data['features'][station_index]['properties']['value']


def get_precip_Locarno_Monti(print_precip=False):
    json_link_precip = 'https://data.geo.admin.ch/ch.meteoschweiz.messwerte-niederschlag-10min/ch.meteoschweiz.messwerte-niederschlag-10min_en.json'
    station_index = 77  # Locarno Monti

    # Get the JSON file
    with urllib.request.urlopen(json_link_precip) as url:
        data = json.loads(url.read().decode())

    # Print, if requested
    if print_precip:
        print('Station: '+data['features'][station_index]['properties']['station_name'])
        print('Time of the measurement: '+data['features'][station_index]['properties']['reference_ts'])
        print('Accumulated precipitation within the last 10 minutes: '+str(data['features'][station_index]['properties']['value'])+' mm')
    
    # Return the accumulated precipitation value
    return data['features'][station_index]['properties']['value']


def get_solar_irrad_Locarno_Monti(print_solar_irr=False):
    json_link_solar_irr = 'https://data.geo.admin.ch/ch.meteoschweiz.messwerte-globalstrahlung-10min/ch.meteoschweiz.messwerte-globalstrahlung-10min_en.json'
    station_index = 78  # Locarno Monti

    # Get the JSON file
    with urllib.request.urlopen(json_link_solar_irr) as url:
        data = json.loads(url.read().decode())

    # Print, if requested
    if print_solar_irr:
        print('Station: '+data['features'][station_index]['properties']['station_name'])
        print('Time of the measurement: '+data['features'][station_index]['properties']['reference_ts'])
        print('Global radiation: '+str(data['features'][station_index]['properties']['value'])+' W/m^2')

    # Return the global radiation value
    return data['features'][station_index]['properties']['value']


def get_time_measurements(print_time=False):
    json_link_T = 'https://data.geo.admin.ch/ch.meteoschweiz.messwerte-lufttemperatur-10min/ch.meteoschweiz.messwerte-lufttemperatur-10min_en.json'
    station_index = 78 # Locarno Monti
    
    # Get the JSON file
    with urllib.request.urlopen(json_link_T) as url:
        data = json.loads(url.read().decode())

    # Print, if requested
    if print_time:
        print('Station: '+data['features'][station_index]['properties']['station_name'])
        print('Time of the measurement: '+data['features'][station_index]['properties']['reference_ts'])

    # Return the time of the measurements
    return data['features'][station_index]['properties']['reference_ts']





if __name__ == '__main__':

    ### Path to the csv file with MeteoSwiss data
    path_csv = 'C:\\xrt\\output\\data\\meteoswiss\\'

    ### Get the latest data from MeteoSwiss
    T = get_T_Locarno_Monti()
    precip = get_precip_Locarno_Monti()
    solar_irr = get_solar_irrad_Locarno_Monti()
    time_meas = get_time_measurements()

    ### Check the existence of the csv file, else create it
    now = datetime.datetime.now()
    csv_filename = 'meteoswiss_' + now.strftime('%Y-%m-%d') + '.csv'
    if not os.path.isfile(path_csv + csv_filename): create_csv_file(path_csv + csv_filename)

    ### Store the MeteoSwiss data in the csv file
    df = pd.read_csv(path_csv + csv_filename)
    df.loc[len(df)] = [time_meas, T, precip, solar_irr]
    df.to_csv(path_csv + csv_filename, index=False)

    ### Add the information to the FIT file
    #hdulist = fits.open(path_folder_fits+file_fit)
    #hdulist[0].header.append(("TEMP-2M", T, "Temperature [degC] 2 meters above ground"))
    #hdulist[0].header.append(("PRECIP", precip, "Accum precip [mm] in the last 10 min"))
    #hdulist[0].header.append(("GLOB-RAD", solar_irr, "Global radiation [W/m^2]"))
    #hdulist[0].header.append(("TIME-ATM", time_meas, "Time of the atmospheric measurements"))

    ### Replace the FIT file with the new one
    #hdulist.writeto(path_folder_fits+file_fit, overwrite=True)


    



