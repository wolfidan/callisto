# -*- coding: utf-8 -*-
"""
Created on 2020-04-13, 2020-05-06
Azimuth and elevation of any fixstar, defined by right-ascension and declination
@author: Christian Monstein
"""
#------------------------------------------------------------------------------

import ephem
import datetime
import numpy as np
import serial

#------------------------------------------------------------------------------

MyHome = ephem.Observer()
MyHome.lon      = '8.7575'    # east +°
MyHome.lat      = '47.205833' # north +°
MyHome.elev     = 414         # altitude in m asl
MyHome.temp     = 20          # °C
MyHome.pressure = 900         # mbar

ra  = 82.5  # Fixstar Taurus A [degree]
dec = 22.01 # Fixstar Taurus A [degree]

MaxRange  = 60.0 # (+/- value) depends on your SAT-rotor type, check data sheet
AziMin    = 120.0
AziMax    = 240.0
EleMin    =  0.0
EleMax    = 90.0

MyComport = 'COM5' # check with device manager (C:\Windows\System32\devmgmt.msc)

#------------------------------------------------------------------------------

#MyHome.date = '2018/06/24 19:11:00' # example for fixed date-time
dt = datetime.datetime.now()
MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
   dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
print ('Current date-time: ',MyHome.date)

fb = ephem.FixedBody();
fb._ra  = np.deg2rad(ra);
fb._dec = np.deg2rad(dec)

fb.compute(MyHome);
alt = np.rad2deg(1.*fb.alt)
azi = np.rad2deg(1.*fb.az);

print("Fixstar AZI= {:4.1f} ELE= {:4.1f}".format(azi,alt))

#------------------------------------------------------------------------------

myazi = -(azi - 180.0) # conversion 0° ... 360° -> +/- 180°
myele = -(alt - 45.0)  # assuming rotor 0° is mounted at 45° elevation

#------------------------------------------------------------------------------

if ((myele>EleMin) and (myele<EleMax) and (myazi>AziMin) and (myazi<AziMax)): # edit for your limitations
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
            
            cmd = 'azi{:6.2f}\r'.format(myazi) 
            DiSEqC.write(cmd.encode()) # set azimuth or hour angle
            
            cmd = 'ele{:6.2f}\r'.format(myele) 
            DiSEqC.write(cmd.encode()) # set elevation or declination
            
            DiSEqC.close()
        
    except IOError:
        print ("Problem communication with tracker. Check COM-port and cables/connectors!")

else:
    print ('Source out of rotor-range of +/-',MaxRange)
    
#------------------------------------------------------------------------------
