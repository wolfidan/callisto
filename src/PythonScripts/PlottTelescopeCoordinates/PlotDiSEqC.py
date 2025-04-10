# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 14:33:50 2016
Generate plots telescope coordinates azimuth, elevationand shows path of geostationary satellites

@author: Chr. Monstein
"""

import numpy as np
from matplotlib import pyplot as plt
import math

#------------------------------------------------------------------------------

path = 'DISEQ-2025-01-07-Sun.txt'

f = open(path, 'r')
sky = np.genfromtxt(f, delimiter=',',skip_header=4)
f.close()

print(sky[0,:])
time = sky[:,0]
azi  = sky[:,1]
ele  = sky[:,2]

#------------------------------------------------------------------------------

def GeoSat(sat_lon):
    lat, lon = 46.337, 8.79 # Locarno Monti Swiss Meteo
    rlat = math.radians(lat)
    rlon = math.radians(lon)
    rsat = math.radians(sat_lon)
    L = rsat - rlon
    D = math.acos(math.cos(rlat) * math.cos(L))
    az = math.degrees(math.acos(-math.tan(rlat) / math.tan(D)))
    az = az if L > 0 else 360 - az
    cd = math.cos(D)
    num = cd - 1 / 6.62
    den = (1 - cd * cd) ** 0.5
    el = math.degrees(math.atan(num / den))
    return(az,el)

#------------------------------------------------------------------------------

plt.figure() # present raw data from spectral verview
plt.plot(azi,ele,'.b',label = 'Telescope -> Sun') 
plt.xlabel("Azimuth [°]")
plt.ylabel("Elevation [°]")
plt.title('Positioning diagram Sun ' + path)
plt.xlim([90,260])
plt.ylim(-12,70)

for sat in range (0,360,5):
    a,e = GeoSat(sat)
    if (e > 0):
        plt.plot(a,e,'.r') 
        plt.plot(360-a,e,'.r') 
plt.plot(a,e,'.r',label='Geo satellites') 

plt.legend(loc="upper left")
plt.savefig('Telescope-Coords.png')

#------------------------------------------------------------------------------

plt.figure() # present raw data from spectral verview
plt.plot(time,azi,'.',label = 'Azimuth',linewidth=1 )
plt.plot(time,ele,'.',label = 'Elevation',linewidth=1) 
plt.xlabel("Time [UT]")
plt.ylabel("[°]")
plt.title('Azimuth/Elevation ' + path)
#plt.xlim([800,1300])
plt.ylim(-12,260)
plt.legend(loc="upper left")
plt.savefig('Telescope-AzEl.png')

#------------------------------------------------------------------------------
