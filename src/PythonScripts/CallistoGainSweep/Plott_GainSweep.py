# -*- coding: utf-8 -*-
"""
Created on 2023-05-18
Data generated with script callisto_gain_sweep.py (twice: onc cold and one hot)
Generate plots hot- and cold measurements with Callisto and estimates Y-factor

@author: Chr. Monstein
"""

import numpy as np
from matplotlib import pyplot as plt

#------------------------------------------------------------------------------

tit = 'WideBand-LNBand-SATF'
f1 = open('GainSweep-2023-11-20-Cold.txt', 'r') # Measurement near zenith (Tcold)
data = np.genfromtxt(f1, delimiter=',',skip_header=5)
f1.close()
pwm   = data[:,0]
Ycold = data[:,1]

f2 = open('GainSweep-2023-11-20-Hot.txt', 'r') # Measurment while looking to ground (Thot = To)
data = np.genfromtxt(f2, delimiter=',',skip_header=5)
f2.close()
Yhot = data[:,1]

#------------------------------------------------------------------------------

plt.figure()# present raw data from spectral verview
plt.plot(pwm,Yhot,'.-',color="red",label="Ground Thot") #o,+,
plt.plot(pwm,Ycold,'.-',color="blue",label="Zenith Tcold") #o,+,
plt.xlabel("Gain control voltage [digit]")
plt.ylabel("Intensity [mV]") # digti is wrong, use mV
plt.title('Raw data from wideband LNB spectral overview')
#plt.text(10870,40,'Receiver band-switching',rotation=90)
plt.legend(loc="lower right")
plt.tick_params(axis='both', which='minor')
#plt.axis([100, 200, 800, 2400])
plt.savefig(tit+'_raw.png',bbox_inches="tight")

#------------------------------------------------------------------------------
plt.figure()
x = pwm
y = (Yhot-Ycold)/25.4 # mV -> dB

poly = np.polyfit(x, y, deg=3) #3...5
plt.plot(x,y,'o', label='data')
plt.plot(x,np.polyval(poly, x), label='fit')
plt.title('Gain sweep for Hot and Cold position to derive pwm-value')
plt.xlabel("Gain control digits [pwm]")
plt.tick_params(axis='both', which='minor')
plt.ylabel("Y-factor [dB]") # digti is wrong, use mV
plt.legend()
plt.show()
plt.savefig(tit+'_YdB.png',bbox_inches="tight")

#------------------------------------------------------------------------------

