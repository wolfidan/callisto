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

path    = 'LC20230904_ADU_NASA.txt'
Thot    = 273.15 + 23.  # 273.15 K + Tamb
Tcold   = 12.  # depends on publication taken into account 12 K ... 37 K
channel = 7  # 1...(7)....10 (prefer a channel without TV-satellites)

with open(path) as file: # read header to get frequency table
    line = file.readline()
header = line.split(',')


f = open(path, 'r') # read data
sky = np.genfromtxt(f, delimiter=',',skip_header=1)
f.close()

time = sky[:,0] # UTC

cali = 1./1024 * 2500./25.4 # 1024 digit/FS, 2500 mV/FS, 25.4 mV/dB
#LCdB = sky[: ,channel]*cali # conversion into dB
#LCdB = sky[: ,1:10]*cali # conversion into dB
LCdB = np.mean(sky[:,1:10], axis=1)/10
LClin = 10.0**(LCdB/10.) # conversion into linear scale


plt.figure()# present raw data from spectral verview
plt.plot(time,LCdB,'-',color='green' )
plt.xlabel("Time [UT]")
plt.ylabel("Intensity [dB]")
plt.title('Light curve Moon transit raw data '+path)
#plt.title('Light curve Sun track&calibrate raw data '+path)
#plt.yscale('log')
plt.ylim(55.6,56)
plt.grid('both', which='minor', color='r', linestyle='-')
plt.grid('both', which='major', color='r', linestyle='-')
plt.tight_layout()
plt.savefig(path+'_'+header[channel]+'_LCmoon.png',bbox_inches='tight')
plt.show()

#------------------------------------------------------------------------------
