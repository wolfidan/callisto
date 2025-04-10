# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 11:09:00 2023

@author: cmons
"""

import numpy as np
import matplotlib.pyplot as plt 

#-------------------------------------------------------------------------------

filename='FluxErrorAll.txt'
header = 3
data   = np.loadtxt(filename,skiprows=header,dtype=str)
datum  = data[:,0]
flux   = data[:,1]
rms    = data[:,2]
lear   = data[:,3]

flux = list(map(float,flux))
rms  = list(map(float,rms))
lear = list(map(float,lear))
print ('Median and rms X-band telescope:   {:6.1f} sfu  {:6.1f} sfu '.format(np.median(flux),np.std(flux)))
print ('Median and rms Learmonth flux:     {:6.1f} sfu  {:6.1f} sfu '.format(np.median(lear),np.std(lear)))
M = np.corrcoef(lear,flux)
print ('Correlation coefficient (Pearson): {:6.3f}'.format(M[0,1]))

fig, ax = plt.subplots(figsize=(12,5))
ax.plot(datum,flux,'-*',label='X-band SRT 07:00-14:00 UT')
ax.plot(datum,lear,'-o',label='Learmonth 05:00 UT')
ax.set_xticks(datum)
ax.set_xticklabels(datum, rotation=90)
ax.set_ylim([300,500])
ax.text(datum[18],326,'[--->',c='blue')
ax.legend(loc="upper left")
ax.set_xlabel('Date')
ax.set_ylabel('SFU')
ax.set_title('Flux comparison X-band SRT vs Learmonth interpolation')
ax.grid()
fig.tight_layout()
fig.savefig("xbandAndLearAll") 
plt.show()

#-------------------------------------------------------------------------------

q1 = 10*np.log10(lear)
q2 = 10*np.log10(flux)
dB = q2 - q1

fig, ax = plt.subplots()
ax.plot(datum,dB,'-x')
ax.set_xticks(datum)

ax.set_xticklabels(datum, rotation=90)
ax.set_ylim([-1.5,1.5])
ax.set_xlabel('Date')
ax.set_ylabel('dB')
ax.set_title('Flux deviation X-band SRT vs Learmonth interpolation')
ax.grid()
fig.tight_layout()
fig.savefig("xbandVsLeardBAll") 
plt.show()

#-------------------------------------------------------------------------------

E = []
N = len(flux)
for i in range(0,N):
    e = (flux[i]-lear[i])/lear[i]*100
    E.append(e)
    
fig, ax = plt.subplots()
ax.plot(datum,E,'-x')
ax.set_xticks(datum)
ax.set_xticklabels(datum, rotation=90)
ax.set_ylim([-30,30])
ax.set_xlabel('Date')
ax.set_ylabel('%')
ax.set_title('Flux deviation X-band SRT vs Learmonth interpolation')
ax.grid()
fig.tight_layout()
fig.savefig("xbandVsLearPercentAll") 
plt.show()

#-------------------------------------------------------------------------------
