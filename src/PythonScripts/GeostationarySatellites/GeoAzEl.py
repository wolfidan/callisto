# -*- coding: utf-8 -*-
"""
Created on 09 Feb 2017
http://www.uhf-satcom.com/uhf/
http://www.rfwireless-world.com/calculators/Antenna-Elevation-Azimuth-Calculator.html
http://www.csgnetwork.com/geosatposcalc.html
@author: Whitham D. Reeve, Anchorage, Alaska USA
Updated: C. Monstein, 2023-03-25

Northern latitudes only!

"""
#------------------------------------------------------------------------------

import numpy as np

#------------------------------------------------------------------------------
# User entries:
latitude  =   46.1723    # North + only
longitude =   8.7875    # West + , East -
satellite =   50.0     # West + , East -

#------------------------------------------------------------------------------

print ('\nTerrestrial coordinate for geostationary satellites')

if satellite < 0:
    satellite = 360.0 + satellite
        
longitude = longitude * np.pi / 180.0
latitude  = latitude  * np.pi / 180.0
satellite = satellite * np.pi / 180.0

delta = (satellite - longitude)

u = np.absolute(delta)

if u == 0:
    azimuth = 0.0
else:
    azimuth = -np.arctan((np.sin(latitude)*(1.0 / np.tan(u)))) + np.pi / 2.0

r = 6373.0 # Earth radius [km]
h = 35752.53 # Geostationary altitude [km]
    
n = np.arccos((np.cos(latitude)) * (np.cos(u)))

elevation = np.arctan((np.cos(n) - r / (r + h)) / np.sin(n))

slant = np.sin(n) * (r + h) / np.cos(elevation)

elevation = elevation * 180.0 / np.pi

azimuth = azimuth * 180.0 / np.pi

if delta > 0.0:
    azimuth = 180.0 + azimuth
else:    
    azimuth = 180.0 - azimuth

el = format(elevation, '.1f')
az = format(azimuth, '.1f')
rg = format(slant, '.0f')

print (' ')
print ('Elevation  : ', el,'°')
print ('Azimuth    : ', az,'°')
print ('Slant Range: ', rg,'km')
    
#---------------------------------------------------------------------------------------