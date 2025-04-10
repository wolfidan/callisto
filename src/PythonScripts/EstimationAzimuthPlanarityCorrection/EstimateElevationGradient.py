# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 14:33:50 2016
Generate correction factor for non-planarity of azimuth plane

@author: Chr. Monstein
"""

import numpy as np
from matplotlib import pyplot as plt

#------------------------------------------------------------------------------

path = 'non-planarity.txt'

f = open(path, 'r')
tilt = np.genfromtxt(f, delimiter=',',skip_header=4)
f.close()

azi = tilt[:,0]

# just for testing
# maybe loop over columns to visualize more than one elevation
#for ii,ele in enumerate(tilt[:,1:]):

ele = tilt[:,2]

#------------------------------------------------------------------------------
EleCoeff = np.polyfit(azi,ele, deg=1) # 1...5
# print ('Coefficients: ',EleCoeff)
x_new = np.linspace(min(azi), max(azi), 100) # 100 = new size
y_fit = np.polyval(EleCoeff, x_new)

st = 'Edit configsun.ini, parameter planecorr = {:6.4f}'.format(EleCoeff[0])
print (st)

#------------------------------------------------------------------------------
ymin = np.min(ele) + 0.06*(np.max(ele-np.min(ele)))
ytext= np.min(ele) + 0.5*(np.max(ele-np.min(ele)))

plt.figure() # present raw data from spectral verview
plt.step(azi,ele,'.-b',label = 'Measured') 
plt.plot(x_new,y_fit,'-r',label = 'Linear fit') 
plt.xlabel("Rotor azimuth [°]")
plt.ylabel("Measured elevation [°]")
plt.title('Error in azimuth plane: ' + path)
plt.text(-40,ymin,st)
plt.text(-75,ytext,'West')
plt.text( 60,ytext,'East')
# plt.xlim([60,260])
# plt.ylim(-12,70)
plt.show()
plt.grid()
plt.legend(loc="upper left")
plt.savefig('AzimuthPlane.png',bbox_inches="tight")

#------------------------------------------------------------------------------
