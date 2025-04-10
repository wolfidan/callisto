# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 19:12:27 2020, 06.05.2020
moves rotors to any position in azimuth and elevation, as long as in range

@author: C. Monstein
"""
import serial


#-----------------------------------------------------------------------------

azi = 200. # target azimuth in degree
ele = 30.  # target elevation in degree

MaxRange  = 60 # (+/- value) depends on your SAT-rotor type, check data sheet
AziMin    = 120.0
AziMax    = 240.0
EleMin    =  0.0
EleMax    = 90.0

MyComport = 'COM5' # check with device manager (C:\Windows\System32\devmgmt.msc)

#------------------------------------------------------------------------------

myazi = -(azi - 180.0) # conversion 0° ... 360° -> +/- 180°
myele = -(ele - 45.0)  # assuming, rotor 0° is mounted at 45° elevation 

if ((ele>EleMin) and (ele<EleMax) and (azi>AziMin) and (azi<AziMax)):
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
    print ('Satellite out of rotor-range of +/-',MaxRange)
    
#------------------------------------------------------------------------------
