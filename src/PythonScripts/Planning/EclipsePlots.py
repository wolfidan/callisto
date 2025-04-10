# -*- coding: utf-8 -*-
"""
Created on 2020-04-13
Updated on 2023-03-25
Azimuth and elevation of sun and moon during eclipse
@author: Christian Monstein
https://nso.edu/eclipse-map-2024/
"""
#------------------------------------------------------------------------------

import ephem
import numpy as np
from matplotlib import pyplot as plt
import math

#------------------------------------------------------------------------------
# Date to be chosen by the user
year  = 2023
month = 10
day   = 14
name  = 'Eclipse Hondo Texas' # title of plot
#------------------------------------------------------------------------------
# Location parameter for selected place

MyHome = ephem.Observer()
MyHome.lon      = '-99.14142'  # east +°
MyHome.lat      = '29.34746'  # north +°
MyHome.elev     = 493         # altitude in m asl
MyHome.temp     = 25          # °C
MyHome.pressure = 900         # mbar

#------------------------------------------------------------------------------

AZI = []
ELE = []
HM  = []
hh  = []

plt.figure()
plt.plot(AZI,ELE)
plt.title("{} on {:02d}-{:02d}-{:04d}".format(name,day,month,year))
plt.xlabel('Azimuth [°]')
plt.ylabel('Elevation [°]')
plt.xticks(np.arange(0, 360+1, 30.0))
plt.xlim(0,360)
plt.ylim(0,90)
plt.text(330,85,'Sun',color='r')
plt.text(330,80,'Moon',color='b')
# for h in range(0,len(hh),3):
#     if ((AZI[h]> 0) and  (AZI[h]<360) and (ELE[h]>0) and (ELE[h]<90)):
#         plt.text(AZI[h],ELE[h],hh[h])
        
sun = ephem.Sun()
moon = ephem.Moon()
dt = []
da = []
de = []
dp = []
for h in range(0,24):
    for m in range(0,60,10):
        #dt = datetime.datetime.now()
   
        MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                year,month,day,h,m,0)
        Tdez = h+m/60
        sun.compute(MyHome)
        moon.compute(MyHome)
        Sazi = math.degrees(sun.az)
        Sele = math.degrees(sun.alt)
        Mazi = math.degrees(moon.az)
        Mele = math.degrees(moon.alt)
        plt.plot(Sazi,Sele,'.r')
        plt.plot(Mazi,Mele,'.b')
        dt.append(Tdez)
        da.append(Sazi-Mazi)
        de.append(Sele-Mele)
        delta = np.sqrt((Sazi-Mazi)**2 + (Sele-Mele)**2)
        dp.append(delta)
        if (delta < 0.03):
            Tazi = Sazi
            Tele = Sele
            Ttim = Tdez
        #print("Sun data: Azimuth   =%6.2f  Elevation   =%6.2f" % (azi-180.0,ele)) 
    if ((h%1)==0):

        if ((Sazi> 0) and  (Sazi<360) and (Sele>0) and (Sele<90)):
            plt.text(Sazi,Sele,h+1) # ?????

plt.grid('both')
plotfilename = '{:4d}-{:02d}-{:02d}_'.format(year,month,day) + name
plt.savefig(plotfilename)
plt.show()

#------------------------------------------------------------------------------

plt.figure()
plt.plot(dt,dp,'.-')
plt.xlim(13,20)
plt.ylim(-1,2)
plt.xlabel('Time [UT]')
plt.ylabel('Angular distance Sun-Moon [°]')
plt.title('Forecast of Totality '+name)
plt.grid('both')
hh = int(Ttim)
mm = int(60*(Ttim - hh))
plt.text(Ttim-1,-0.15,"Totality: ~{:02d}:{:02d} UT".format(hh,mm))
plt.text(Ttim-1,-0.35,"Azimuth: {:04.1f}°".format(Tazi))
plt.text(Ttim-1,-0.55,"Elevation: {:04.1f}°".format(Tele))
plt.savefig('Totality.png')
plt.show()

#------------------------------------------------------------------------------
