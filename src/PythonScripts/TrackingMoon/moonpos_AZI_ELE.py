# -*- coding: utf-8 -*-
"""
This script tracks the Moon in azimuthal mode (azimuth and elevation)
Carefully check paramter in MyLocation
It is sufficient to execute this script once every minute
Created on Thu Sep 14 20:13:50 2017
Updated: 09.05.2020, 05.05.2023, 13.05.2023, 23.06.2023
@author: Chr. Monstein
"""
# http://rhodesmill.org/pyephem/tutorial.html
# http://rhodesmill.org/pyephem/quick#other-observer-methods
#------------------------------------------------------------------------------

import ephem
import serial
import datetime
import math
import numpy as np
import os.path
import configparser

#------------------------------------------------------------------------------
# Read configuration file
config = configparser.ConfigParser()
config.read('configmoon.ini')
  
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
EleGnd   = config.getfloat('Calibration','EleGnd')
dEleCali = config.getfloat('Calibration','dEleCali')

#------------------------------------------------------------------------------
# Scanning parameter
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]

Hot  = np.array(getlist(config.get('TelescopeSwitching','Hot' )), dtype=np.int)
Cold = np.array(getlist(config.get('TelescopeSwitching','Cold')), dtype=np.int)
Src  = np.array(getlist(config.get('TelescopeSwitching','Src' )), dtype=np.int)

#------------------------------------------------------------------------------

TelescopeParked = False # to prevent endless parking in the night

#------------------------------------------------------------------------------

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
            
            cmd = 'azi{:6.2f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            
            cmd = 'ele{:6.2f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            
            DiSEqC.close()
        
    except IOError:
        print ("Problem communication with tracker!")
    
#------------------------------------------------------------------------------

def Log(x,y,msg):
    dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
    filename ='DISEQ-'+dat+'-Moon.txt'
    if (os.path.exists(filename)):
        with open(filename, "a") as fp:   # Save x/y-data in file
            hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            # add decimal time
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
            #st = hms +',{:8.3f},'.format(x) + '{:8.3f}, '.format(y) + msg
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
            #st = hms +',{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
            fp.write(st + '\n') 

#------------------------------------------------------------------------------

dt = datetime.datetime.now() # PC must run on UT or GPS-time
MyLocation.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
   dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
#MyLocation.date = '2023/03/14 07:24:00' # example for testing at high noon
print ('Selected date-time: ',MyLocation.date,' UT')

moon = ephem.Moon()
moon.compute(MyLocation)
azi  = math.degrees(moon.az)
ele0 = math.degrees(moon.alt)
print("Moon data: Azimuth   = %7.2f  Elevation   = %6.2f" % (azi,ele0))

lst = MyLocation.sidereal_time()
ha = (lst - moon.ra)/math.pi*180.0
dec = math.degrees(moon.dec)
print("Moon data: Hourangle = %7.2f  Declination = %6.2f" % (ha,dec))

# Telescope switching mode evaluation
if (dt.minute in Hot):
    ele = EleGnd
    Log(azi,ele,'Hot reference on ground (To)')
else:
    if (dt.minute in Cold):
        ele = ele0 + dEleCali
        Log(azi,ele,'Off source in elevation')
    else:
        if (dt.minute in Src):
            ele = ele0
            Log(azi,ele,'On source in elevation')

#------------------------------------------------------------------------------
# Correction of non-planarity

def ElevationRevision(myazi,myele):
    # constant value depends on local conditions and may change
    # constant can be found by measureing elevation at all azimuth by linear regression
    return (myele - 0.0184*myazi)

#------------------------------------------------------------------------------

# tracking block

RotorAzi = azidir*(azi - aziref) # conversion 0° ... 360° -> +/- 180°
RotorEle = eledir*(ele - eleref)
RotorEle = ElevationRevision(RotorAzi,RotorEle) # +/-1.25°

if ((azi > AziMax) or (azi < AziMin) or 
    (ele < EleMin) or (ele > EleMax) or
    (np.abs(RotorEle) > MaxRange)    or
    (np.abs(RotorAzi) > MaxRange)):
    if (TelescopeParked == False):
        Park() # no more tracking activity useful
        print("Rotor position or Moon out of range")
        Log(RotorAzi,RotorEle,"Rotor position possibly out of range")
        Log(azi,ele,"Moon position possibly out of range")
        TelescopeParked = True
else:
    print("Rotor and Moon position in range: Azimuth   = %7.2f  Elevation   = %6.2f" % (RotorAzi,RotorEle))
    Log(azi,ele,'New tracking position')   
            
        
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
            
            cmd = 'ele{:6.2f}\r'.format(RotorEle) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            
            cmd = 'azi{:6.2f}\r'.format(RotorAzi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            
            DiSEqC.close()
        
    except IOError:
        print ("Problem communication with elevation tracker!")

#------------------------------------------------------------------------------