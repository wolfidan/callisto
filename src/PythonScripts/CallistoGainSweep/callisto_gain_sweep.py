# -*- coding: utf-8 -*-
"""
Created on Fri Dec 08 11:30:50 2017
Sweeps through all gain values of Callisto, once for hot = ground and once for cold = sky
to determine best gain value whihc needs to be entered in file callisto.cfg
Goal: maximum sensitivity = lowest noise figure but, reduce risk of saturation

@author: Monstein
"""
#---------------------------------------------------------------------------

import time
import serial
from matplotlib import pyplot as plt
import re
import datetime
import os.path

#---------------------------------------------------------------------------

AGC_start =  0 # lower end of the spectrum [MHz]
AGC_stop  =  250 # high end of the spectrum [MHz]
AGC_step  =   5  # frequency step [MHz]
frequency = 10931.0 # LNB frequency MHz
LO = 10410.0 # local oscillator MHz
debug   =  False  # enable or disable text printing, default True
title = 'Hot'

#---------------------------------------------------------------------------

def Log(x,y):
    dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
    filename ='GainSweep-'+dat+'-'+title+'.txt'
    if (os.path.exists(filename)):
        with open(filename, "a") as fp:   # Save x/y-data in file
            #hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
            # add decimal time
            st = '{:03.0f},'.format(x) + '{:03.0f}'.format(y) 
            fp.write(st+'\n') 
    else:
        with open(filename, "a") as fp:   # Save x/y-data in file
            fp.write('Callisto gain sweep\n') # write header information
            fp.write('Author: Christian Monstein, HB9SCT\n') # write header information
            fp.write('Version: 2023-05-18\n') # write header information
            fp.write('\n')
            fp.write('pwm[digit], Voltage[mV]\n') # write header information

#---------------------------------------------------------------------------
# Initialization serial port
try:
    callisto = serial.Serial(
         port     = 'COM7',
         baudrate = 115200,
         bytesize = serial.EIGHTBITS,
         parity   = serial.PARITY_NONE,
         xonxoff  = False,
         rtscts   = False,
         dsrdtr   = False, 
         timeout  = 0)

    callisto.isOpen()    
    print("Connected to Callisto: " + callisto.portstr)
except IOError:
    callisto.close()
    callisto.open()
    print("Callisto was already open, was closed and opened again!")
#---------------------------------------------------------------------------
# Transmit paramter format and pwm (gain)   
cmd = 'F0{:07.1f}\r'.format(frequency-LO)
callisto.write(cmd.encode())

cmd='%5\r' # formatting output
callisto.write(cmd.encode()) # formatting output

cmd='O0\r' # pwm
callisto.write(cmd.encode()) # pwm

time.sleep(2.5) # wait until gain is stable

#---------------------------------------------------------------------------
# empty potential garbage in the serial buffer
for c in callisto.read():
    #print (c)
    if c == '':
        break
#---------------------------------------------------------------------------
# # Define frequency list    
# def my_range(AGC_start, AGC_stop, AGC_step):
#     while AGC_start <= AGC_stop:
#         yield AGC_start
#         AGC_start += AGC_step
# #---------------------------------------------------------------------------
# Trigger spectrometer and await for detector voltage [mV]
spec_v  = []
spec_x  = []
for agc in range(AGC_start,AGC_stop+AGC_step,AGC_step):
    cmd = 'O{:03.0f}\r'.format(agc) 
    callisto.write(cmd.encode())
    time.sleep(0.5) # charge time
    cmd = 'A0\r'
    callisto.write(cmd.encode())
    time.sleep(0.1) # allow time to measure and transmit     
    out = callisto.readline()
    text = str(out, 'utf-8')
    if (("$CRX:" in text) and ((len(text)==18) or (len(text)==17))):
        if (debug):
            print (text)
        r = re.compile('[:,]')
        vv = r.split(text)
        volt = float(vv[2])
        spec_v.append(volt)
        spec_x.append(agc)
        Log(agc,volt)
        print (agc,volt)
        
callisto.close()

#---------------------------------------------------------------------------
# Plot raw data
plt.figure(figsize=(8,5)) # calculate external rfi and plot it
plt.grid()
plt.xlabel("Gain control [digit]",fontsize=15)
plt.ylabel("Intensity [mV]",fontsize=15)
plt.title('Callisto gain sweep '+title,fontsize=15)
#plt.ylim([500,2400])
plt.tick_params(labelsize=14)
plt.plot(spec_x,spec_v,'-r.')
plt.savefig('Callisto_mV_'+title+'.png')

#---------------------------------------------------------------------------
