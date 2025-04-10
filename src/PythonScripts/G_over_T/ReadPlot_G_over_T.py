# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 14:33:50 2016
Updated on July 7 19:22
Generate plots of external rfi from spectral overview files OVS

@author: Chr. Monstein
"""

import numpy as np
from matplotlib import pyplot as plt

path = 'C:/CALLISTO-01/OVSfiles/'
LO = 10.41 # GHz

f = open(path+'OVS_SWISS-METEO_X-Band_20231120_132450_01.prn', 'r')
Dcold = np.genfromtxt(f, delimiter=';',skip_header=1)
f.close()

f = open(path+'OVS_SWISS-METEO_X-Band_20231120_132019_01.prn', 'r')
Dhot = np.genfromtxt(f, delimiter=';',skip_header=1)
f.close()
freq = Dcold[:,0]/1000. + LO
Icold = Dcold[:,1]
Ihot  = Dhot [:,1]
ydB = (Ihot - Icold)/25.4

Tcold = 12 # from literature
Thot = 273.15 + 10.0 # from outdoor sensor
Y = 10**(np.median(ydB)/10)
Trx = (Thot - Y*Tcold) / (Y - 1)
GdB = 36
G = 10**(GdB/10)
GoverT_dB = 10*np.log10(G/Trx)
print(G)
print(Y)
print(Trx)
print('G/T: %s dB' % GoverT_dB)

plt.figure()# present raw data from spectral verview
plt.plot(freq,Ihot,'-',color="red",label='Hot ground') #o,+,
plt.plot(freq,Icold,'-',color="blue",label='Cold sky') #o,+,
plt.xlabel("Frequency [GHz]")
plt.ylabel("Intensity [mV]")
plt.title('Raw data from spectral overview ')
plt.legend(loc="upper left")
plt.savefig(path+'Xband_raw.png',bbox_inches='tight')


plt.figure() # calculate external rfi and plot it
#plt.figure(figsize=(20,10)) # calculate external rfi and plot it
plt.plot(freq,ydB,'-', color="red")
plt.grid()
plt.xlabel("Frequency [GHz]")
plt.ylabel("Y-factor [dB]")
plt.title('System response hot ground - Sky')
plt.savefig(path+'Xband-Y.png',bbox_inches='tight')

