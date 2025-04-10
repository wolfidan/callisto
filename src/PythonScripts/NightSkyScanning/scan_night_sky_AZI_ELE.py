# -*- coding: utf-8 -*-
"""
This script (strongly based on sunpos_AZI_ELE.py) scans the night sky in 
azimuthal mode (azimuth and elevation). Carefully check parameter in external file 
config.ini. It is sufficient to execute this script once every minute, during the 
scan interval. Check the array of scans called coordinates_night-scan-plan.csv

Created on Tue Mar 04
@author: Andrea Francesco Battaglia
"""

#-----------------------------------------------------

import ephem
import serial # was serial
import datetime
#import math
import os.path
import numpy as np
import time
import configparser
import pandas as pd

#-----------------------------------------------------

# Absolute path to the csv file containing the times and coordinates of the scan
csv_coordinates_scan1 = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_TV-sat.csv'
csv_coordinates_scan2 = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_zenith-ground.csv'
csv_coordinates_scan3 = 'C:\\xrt\\src\\PythonScripts\\NightSkyScanning\\ScanningModes\\coordinates_night-scan_full-sky.csv'

#-----------------------------------------------------

# Read configuration file
config = configparser.ConfigParser()
config.read('C:\\xrt\\src\\PythonScripts\\TrackingSun\\configsun.ini')

# Location parameter
MyLocation          = ephem.Observer()
MyLocation.lon      = config.get('Location','longitude')
MyLocation.lat      = config.get('Location','latitude')
MyLocation.elev     = config.getfloat('Location','elevation')
MyLocation.temp     = config.getfloat('Location','temperature')
MyLocation.pressure = config.getfloat('Location','pressure')

# Tracker parameter
MaxRange  = config.getfloat('Tracker','MaxRange')
AziMin    = config.getfloat('Tracker','AziMin')
AziMax    = config.getfloat('Tracker','AziMax')
EleMin    = config.getfloat('Tracker','EleMin')
EleMax    = config.getfloat('Tracker','EleMax')
AziPark   = config.getfloat('Tracker','AziPark')
ElePark   = config.getfloat('Tracker','ElePark')
aziref    = config.getfloat('Tracker','aziref')
eleref    = config.getfloat('Tracker','eleref')
azidir    = config.getfloat('Tracker','azidir')
eledir    = config.getfloat('Tracker','eledir')
planecorr = config.getfloat('Tracker','planecorr')
MyComport = config.get('Tracker','MyComport')

# Calibration parameter
#CaliProc = config.getboolean('Calibration','CaliProc')
#AziCold  = config.getfloat('Calibration','AziCold')
#EleCold  = config.getfloat('Calibration','EleCold')
#AziHot   = config.getfloat('Calibration','AziHot')
#EleHot   = config.getfloat('Calibration','EleHot')
#TimeCali = config.getint('Calibration','TimeCali')

#-----------------------------------------------------

# Scanning parameter
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]

scanning = config.getboolean('Scanning','scanning')
Dazi = np.array(getlist(config.get('Scanning','Dazi')), dtype=np.float32)
Dele = np.array(getlist(config.get('Scanning','Dele')), dtype=np.float32)

#-----------------------------------------------------

TelescopeParked = False # to prevent endless parking during the night

#-----------------------------------------------------

# Park telescope after sun-set or after source is out of view
def Park():
    RotorAzi = azidir*(AziPark - aziref) # conversion 0° ... 360° -> +/- 180°
    RotorEle = eledir*(ElePark - eleref)
    print ("Parking at Azimuth  = %7.2f  Elevation   = %6.2f" % (AziPark,ElePark))
    
    try:
        DiSEqC = serial.Serial(
             port     = MyComport,
             baudrate = 9600,
             bytesize = serial.EIGHTBITS,
             parity   = serial.PARITY_NONE,
             timeout  = 2)
        
        print(DiSEqC)
             
        if (DiSEqC.isOpen()):
            print ("Successfully connected to antenna tracker at: "+DiSEqC.portstr)
            cmd = 'max{:6.2f}\r'.format(MaxRange) 
            DiSEqC.write(cmd.encode()) # set MaxRange
            time.sleep(0.02)

            cmd = 'azi{:8.4f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            time.sleep(0.02)

            cmd = 'ele{:8.4f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            time.sleep(0.02)
            
            DiSEqC.close()
        
    except Exception as ee: #IOError:
        print("Problem communication with tracker!")
        print(ee)
        raise

#-----------------------------------------------------

# Logging events
def Log(x,y,msg):
    dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
    filename ='DISEQ-' + dat + '_nigh-sky-scan.txt'
    if (os.path.exists(filename)):
        with open(filename, "a") as fp:   # Save x/y-data in file
            hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(x) + '{:8.3f}, '.format(y) + msg
            fp.write(st+'\n') 
    else:
        with open(filename, "a") as fp:   # Save x/y-data in file
            fp.write('DiSEqC radio telescope tracker\n') # write header information
            fp.write('Author: Christian Monstein, HB9SCT\n') # write header information
            fp.write(' - Modified by Andrea F. Battaglia for the night-sky scan\n') # write header information
            fp.write('Version: 2025-03-03\n') # write header information
            fp.write('\n')
            hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
            fp.write('TIME [UT], Azimuth[deg], Elevation[deg], Message\n') # write header information
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
            fp.write(st) 

#-----------------------------------------------------

# Correction of non-planarity
def ElevationRevision(myazi,myele):
    # constant value depends on local conditions and may change
    # constant can be found by measureing elevation at all azimuth by linear regression
    return (myele - planecorr*myazi)

#-----------------------------------------------------

# Ephemerides evaluiation
dt = datetime.datetime.now() # PC must run on UT or GPS-time
MyLocation.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
   dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
#MyLocation.date = '2023/03/14 07:24:00' # example for testing at any date&time
print ('Selected date-time: ',MyLocation.date,' UT')

### Since the scanning goes over two days, and we need a criterion for the
## parking of the telescope, we had to do this trick (not so elegant...)
now = datetime.datetime.now()
print("Current date and time: ", now)
if now.hour >= 19 and now.hour < 24:
    date_start = now.strftime("%Y-%m-%d")
    date_end = now + datetime.timedelta(days=1)
    date_end = date_end.strftime("%Y-%m-%d")
    print("Start date: ", date_start)
    print("End date: ", date_end)
elif now.hour >= 0 and now.hour < 7:
    date_end = now.strftime("%Y-%m-%d")
    date_start = now - datetime.timedelta(days=1)
    date_start = date_start.strftime("%Y-%m-%d")
    print("Start date: ", date_start)
    print("End date: ", date_end)

### Get the times and coordinates of the scan
df_scan1 = pd.read_csv(csv_coordinates_scan1)
df_scan1['time_UT'] = date_start + ' ' + df_scan1['time_UT']
df_scan2 = pd.read_csv(csv_coordinates_scan2)
df_scan2['time_UT'] = date_start + ' ' + df_scan2['time_UT']
df_scan3 = pd.read_csv(csv_coordinates_scan3)
df_scan3['time_UT'] = date_end + ' ' + df_scan3['time_UT']
df_scan = pd.concat([df_scan1, df_scan2, df_scan3], ignore_index=True)
times_scan = df_scan['time_UT'].values  
azi_scan = df_scan['azimuth_deg'].values
ele_scan = df_scan['elevation_deg'].values

### Get closes index to now
idx_now = np.argmin(np.abs(pd.to_datetime(times_scan) - pd.to_datetime(now)))
times_scan[idx_now]
azi = azi_scan[idx_now]
ele = ele_scan[idx_now]

#-----------------------------------------------------

# Scan the night sky
RotorAzi = azidir*(azi - aziref) # conversion 0° ... 360° -> +/- 180°
RotorEle = eledir*(ele - eleref) # conversion 0° ... 90° -> +/- 45°
RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°

### Since the coordinates should be already in the correct range,
## as a criterion to stop the scanning is the time.
if pd.to_datetime(now) > (pd.to_datetime(times_scan[-1]) + pd.Timedelta('1 minute')):
    if (TelescopeParked == False):
        #Park() # no more tracking activity useful
        print('Night-sky scan ended')
        Log(0,0,'Night-sky scan ended')
        TelescopeParked = True
else:
    print("Rotor and sun position in range: Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))
    
    try:
        DiSEqC = serial.Serial(
             port     = MyComport,
             baudrate = 9600,
             bytesize = serial.EIGHTBITS,
             parity   = serial.PARITY_NONE,
             timeout  = 2)
        
        print(DiSEqC)
             
        if (DiSEqC.isOpen()):
            print ("Successfully connected to antenna tracker at: "+DiSEqC.portstr)
            cmd = 'max{:6.2f}\r'.format(MaxRange) 
            DiSEqC.write(cmd.encode()) # set MaxRange
            time.sleep(0.02)
            cmd = 'ele{:8.3f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            time.sleep(0.02)

            cmd = 'azi{:8.3f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            time.sleep(0.02)

            Log(azi,ele,'New tracking position')   

            ### The following is no more needed
            '''
            if (CalibrationActive == False):
                if (scanning == True):
                    print ('Start scanning')
                    for i in range(0,len(Dazi)):
                        RotorAzi = azidir*(azi - aziref + Dazi[i]) # conversion 0° ... 360° -> +/- 180°
                        cmd = 'azi{:8.3f}\r'.format(RotorAzi) 
                        DiSEqC.write(cmd.encode()) # set azimuth or hour angle
                        time.sleep(1.1)
                        
                        RotorEle = eledir*(ele - eleref + Dele[i])
                        RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°
                        cmd = 'ele{:8.3f}\r'.format(RotorEle) 
                        DiSEqC.write(cmd.encode()) # set elevation or declination
                        time.sleep(2.2)
                        print("New scanning position Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))                    
                        Log(azi + Dazi[i],ele + Dele[i],'New scanning position')   
            '''
            DiSEqC.close()
        
    except IOError:
        print ("Problem communication with elevation tracker!")
        raise

#-----------------------------------------------------