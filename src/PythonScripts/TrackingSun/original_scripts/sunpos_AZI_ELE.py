# -*- coding: utf-8 -*-
"""
This script tracks the Sun in azimuthal mode (azimuth and elevation)
Carefully check parameter in external file config.ini
It is sufficient to execute this script once every minute
Created on Thu Sep 14 20:13:50 2017
Updated: 09.05.2020, 05.05.2023, 13.05.2023, 04.08.2023, 23.08.2023
@author: Chr. Monstein
"""
#------------------------------------------------------------------------------

import ephem
import serial
import datetime
import math
import os.path
import numpy as np
import time
import configparser

#------------------------------------------------------------------------------
# Read configuration file
config = configparser.ConfigParser()
config.read('configsun.ini')

# with open(os.path.join('C:\\xrt\\src\\PythonScripts\\TrackingSun','test.txt'),'w') as f:
#     f.write('hello!')
  
# Location parameter
MyLocation          = ephem.Observer()
MyLocation.lon      = config.get('Location','longitude')
MyLocation.lat      = config.get('Location','latitude')
MyLocation.elev     = config.getfloat('Location','elevation')
MyLocation.temp     = config.getfloat('Location','temperature')
MyLocation.pressure = config.getfloat('Location','pressure')

#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Calibration parameter
CaliProc = config.getboolean('Calibration','CaliProc')
EleGnd   = config.getfloat('Calibration','EleGnd')
dEleCali = config.getfloat('Calibration','dEleCali')
TimeCali = config.getint('Calibration','TimeCali')

#------------------------------------------------------------------------------
# Scanning parameter
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]

scanning = config.getboolean('Scanning','scanning')
Dazi = np.array(getlist(config.get('Scanning','Dazi')), dtype=np.float32)
Dele = np.array(getlist(config.get('Scanning','Dele')), dtype=np.float32)

#------------------------------------------------------------------------------

TelescopeParked = False # to prevent endless parking during the night

#------------------------------------------------------------------------------
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
        print(ee)
        print ("Problem communication with tracker!")
    
#------------------------------------------------------------------------------
# Logging events

def Log(x,y,msg):
    dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
    filename ='DISEQ-' + dat + '-Sun.txt'
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
            fp.write('Version: 2023-04-27\n') # write header information
            fp.write('\n')
            hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
            fp.write('TIME [UT], Azimuth[°], Elevation[°], Message\n') # write header information
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
            fp.write(st + '\n') 

#------------------------------------------------------------------------------
# Ephemerides evaluiation

dt = datetime.datetime.now() # PC must run on UT or GPS-time
MyLocation.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
   dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
#MyLocation.date = '2023/03/14 07:24:00' # example for testing at any date&time
print ('Selected date-time: ',MyLocation.date,' UT')

sun = ephem.Sun()
sun.compute(MyLocation)
azi = math.degrees(sun.az)
ele = math.degrees(sun.alt)
print("Sun data: Azimuth   = %7.2f  Elevation   = %6.2f" % (azi,ele))

lst = MyLocation.sidereal_time()
ha = (lst - sun.ra)/math.pi*180.0
dec = math.degrees(sun.dec)
print("Sun data: Hourangle = %7.2f  Declination = %6.2f" % (ha,dec))

CalibrationActive = False
if (CaliProc): # find out which position for calibration
    if ((dt.minute%TimeCali) == 0): # calibration cold sky up or down
        ele = ele + dEleCali
        Log(azi,ele,'Cold sky reference deviation in elevation')
        CalibrationActive = True
        
    if (((dt.minute-1)%TimeCali) == 0): # calibration hot ground at low elevation
        ele = EleGnd
        Log(azi,ele,'Hot reference on ground (To)')
        CalibrationActive = True
        
    if (((dt.minute-2)%TimeCali) == 0): # Calibration finished, back to the sun
        Log(azi,ele,'Back to the Sun') # keep original position of sun
        CalibrationActive = True

#------------------------------------------------------------------------------
# Correction of non-planarity

def ElevationRevision(myazi,myele):
    # constant value depends on local conditions and may change
    # constant can be found by measureing elevation at all azimuth by linear regression
    return (myele - planecorr*myazi)

#------------------------------------------------------------------------------
# track-/scan-/calibration-block

RotorAzi = azidir*(azi - aziref) # conversion 0° ... 360° -> +/- 180°
RotorEle = eledir*(ele - eleref) # conversion 0° ... 90° -> +/- 45°
RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°

if ((azi > AziMax) or (azi < AziMin) or 
    (ele < EleMin) or (ele > EleMax) or
    (np.abs(RotorEle) > MaxRange)    or
    (np.abs(RotorAzi) > MaxRange)):
    if (TelescopeParked == False):
        Park() # no more tracking activity useful
        print("Rotor position or sun out of range")
        Log(RotorAzi,RotorEle,"Rotor position possibly out of range")
        Log(azi,ele,"Sun position possibly out of range")
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
            DiSEqC.close()
        
    except Exception as ee: #IOError:
        print(ee)
        print ("Problem communication with elevation tracker!")

#------------------------------------------------------------------------------