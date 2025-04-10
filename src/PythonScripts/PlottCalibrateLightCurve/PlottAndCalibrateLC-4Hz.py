# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 14:33:50 2016
Generate plots of Callisto light curve(s) and calibrate
Works only, when FIT-file generation starts at T%15 minutes
Requires integration factor inside frq?????.cfg of 40 = 40*0.25s = 10 seconds time resolution in LC
@author: C. Monstein
"""

import numpy as np
from matplotlib import pyplot as plt

#------------------------------------------------------------------------------

path    = 'LC20231121_ADU_X-Band.txt'
Thot    = 273.15 + 10.  # 273.15 K + Tamb
Tcold   = 12.  # depends on publication taken into account 12 K ... 37 K
channel = 7  # 1...(7)....10 (prefer a channel without TV-satellites)
Tstart  = 9.50 # start calibration > sun-rise, used to mask data
Tstop   = 9.75 # stop calibration <  sun-set, used to mask data
dT      = 0.25 # calibration period in hours = length of a FIT-file
LearmonthFlux = 358.7 # check here: https://owenduffy.net/calc/qsrf/qsrfr.php


with open(path) as file: # read header to get frequency table
    line = file.readline()
header = line.split(',')


f = open(path, 'r') # read data
sky = np.genfromtxt(f, delimiter=',',skip_header=1)
f.close()

time = sky[:,0] # UTC

cali = 1./1024 * 2500./25.4 # 1024 digit/FS, 2500 mV/FS, 25.4 mV/dB
LCdB = sky[: ,channel]*cali # conversion into dB
LClin = 10.0**(LCdB/10.) # conversion into linear scale

N = int((Tstop-Tstart)/dT) # number of calibration steps
print ('Number of analysis steps ',N)
 

plt.figure()# present raw data from spectral verview
plt.plot(time,LClin,'-',color='green',label = header[channel] )
plt.xlabel("Time [UT]")
plt.ylabel("Intensity [linear]")
plt.title('Light curve raw data '+path)
plt.yscale('log')
plt.grid('both', which='minor', color='r', linestyle='-')
plt.grid('both', which='major', color='r', linestyle='-')
plt.legend(loc="upper left")
plt.savefig(path+'_'+header[channel]+'_LC.png',bbox_inches='tight')
plt.show()

#------------------------------------------------------------------------------

mask = (time > Tstart) * (time < Tstop) # limit analysis to the range of no obstructions

time = time[mask]
LClin = LClin[mask]
M = len(LClin) # number of measurement points within selected time-range

Acold = []
Ahot  = []
j = 0
i = 0
StartCold = 15.0 # seconds after calibration command. Find out by looking into details of light-curve
StopCold  = 57.0 # seconds after calibration command. Find out by looking into details of light-curve
while (i < M):
    t = time[i]
    y = LClin[i]
    T1 = Tstart + j*dT + StartCold / 3600 # 22 start cold sky
    T2 = Tstart + j*dT + StopCold  / 3600 # 51 stop cold sky
    Icold = 0
    k = 0
    if ((t >= T1) and (t <= T2)):
        Icold = Icold + y
        k = k + 1
        print ('i, t, y: ',i,t,y)

        Icold = Icold / k
        print ('Icold, k: ',Icold,k,'\n')
        Acold.append(Icold)
        j = j + 1
    else:
        i = i + 1
    
j = 0
i = 0
StartHot = 85. + StartCold # seconds after calibration command. Find out by looking into details of light-curve
StopHot  = 115. + StopCold  # seconds after calibration command. Find out by looking into details of light-curve
while (i < M):
    t = time[i]
    y = LClin[i]
    T1 = Tstart + j*dT + StartHot / 3600 # 82 start hot ground
    T2 = Tstart + j*dT + StopHot  / 3600 # 125 stop hot ground
    Ihot = 0
    k = 0
    if ((t >= T1) and (t <= T2)):
        Ihot = Ihot + y
        k = k + 1
        print ('i, t, y: ',i,t,y)

        Ihot = Ihot / k
        print ('Ihot, k: ',Ihot,k,'\n')
        Ahot.append(Ihot)
        j = j + 1
    else:
        i = i + 1

#------------------------------------------------------------------------------

plt.figure()
plt.plot(time,LClin, color='black',label = header[channel])
plt.xlabel("Time [UT]")
plt.ylabel("Intensity [linear]")
plt.title('Hot and cold measurement '+path)
plt.yscale('log')
plt.legend(loc="lower right")
plt.savefig(path+'_'+header[channel]+'_Interpolation.png',bbox_inches='tight')
plt.show()

x=np.arange(0,len(LClin))

plt.figure()
plt.plot(x,LClin, color='black',label = header[channel])
plt.xlabel("Time [pixel number]")
plt.ylabel("Intensity [linear]")
plt.title('Hot and cold measurement '+path)
plt.yscale('log')
plt.legend(loc="lower right")
plt.savefig(path+'_'+header[channel]+'_Interpolation-Pixel.png',bbox_inches='tight')
plt.show()

#------------------------------------------------------------------------------
T5 = int(15) # start Isun scan, pixel
T6 = int(90)
sun = LClin[T5:T6]
sunpeak = np.max(sun)
print ('Peak sun {:8.1f} digit '.format(sunpeak))

i = 0
I = 0
for isun in sun:
    if (isun> 0.97*sunpeak):
        I = I + isun
        i = i + 1
Isun = I/i
print ('Average sun (0.97*SunPeak) {:8.1f} digit from {:3.0f} measured pixels'.format(Isun,i))

IhotX  = np.mean(Ahot)
IcoldX = np.mean(Acold)

Tsun  = (Thot-Tcold)/(IhotX-IcoldX) * (Isun-IcoldX) + Tcold
print ('Antenna temperature: {:5.1f} K'.format(np.median(Tsun))) # observed flux

Tgas = Tsun*(3.3/0.5)**2
print ('Corona gas temperature due to beam dilution: {:6.0f} K'.format(np.median(Tgas))) # observed flux

#------------------------------------------------------------------------------

freqstring = header[channel] # selected frequency channel
frequency  = float(freqstring[0:9]) # Frequency in MHz
Lambda = 3e8/(frequency * 1e6) # Wavelength in meter
Bk     = 1.38e-23 # Boltzmann constant in J/K
Gain   = 10**3.6 # 36 dB antenna gain
sfu    = 1e22 # conversion flux into sfu
Ssun   = 8*np.pi*Bk*Tsun / (Lambda**2 * Gain) * sfu
Smed = np.median(Ssun)
print ('Frequency: ' + freqstring)
print ('Median flux:                           {:6.1f} sfu'.format(Smed)) # observed flux
print ('Learmonth interpolation of:            {:6.1f} sfu'.format(LearmonthFlux)) # reference flux

Sback = (2*Bk*Tcold)/(Lambda**2*Gain)*4*np.pi * sfu
print ('Background flux:                       {:6.1f} SFU'.format(Sback))
print ('Final flux above background:           {:6.1f} SFU'.format(Smed-Sback))

Squiet = 2.79e-5 * frequency**1.748 # quiet Sun flux after A. O. Benz
print ('Benz quiet flux:                       {:6.1f} sfu'.format(Squiet))

# #------------------------------------------------------------------------------
