# -*- coding: utf-8 -*-
"""
Created on 2020-04-13, 2020-05-06, 2023-11-09
Azimuth and elevation of any fixstar, defined by right-ascension and declination at given date
Azimuth and elevation sun and path of geostationary satellites
@author: Christian Monstein
"""
#------------------------------------------------------------------------------

import ephem
import numpy as np
from matplotlib import pyplot as plt
import math

#------------------------------------------------------------------------------
# Date to be chosen by the user
year  = 2024
month =  12
day   = 17
#------------------------------------------------------------------------------
# Location parameter for Locarno Monti

MyHome = ephem.Observer()
MyHome.lon      = '8.7875'    # east +°
MyHome.lat      = '46.1723'  # north +°
MyHome.elev     = 374       # altitude in m asl
MyHome.temp     = 15        # °C
MyHome.pressure = 900       # mbar

#------------------------------------------------------------------------------
# Telescope parameter for dishes 5m and 7m

AziMin    =  180-77 # east rotor limit
AziMax    =  180+77 # west rotor limit
EleMin    =  0 # 5m safety
EleMax    =  85.0
xT = [AziMin,AziMax,AziMax,AziMin,AziMin] # telescope frame azimuth
yT = [EleMin,EleMin,EleMax,EleMax,EleMin] # telescope frame elevation

#------------------------------------------------------------------------------

def GeoSat(sat_lon):
    lat, lon = math.degrees(float(MyHome.lat)), math.degrees(float(MyHome.lon))
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

def MyAstro(ra,dec,name):
    fb = ephem.FixedBody();
    fb._ra  = np.deg2rad(ra);
    fb._dec = np.deg2rad(dec)
    
    AZI = []
    ELE = []
    HM  = []
    hh  = []

    for h in range(0,24):
        for m in range(0,60,20):
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year,month,day,h,m,0)
            
            fb.compute(MyHome);
            alt = np.rad2deg(1.*fb.alt)
            azi = np.rad2deg(1.*fb.az);
            hm = h+m/60
            HM.append(hm)
            AZI.append(azi)
            ELE.append(alt)
            st = '{:4.0f}'.format(hm)
            hh.append(st)

    
    plt.figure()
    plt.plot(AZI,ELE,'.b')
    plt.plot(xT,yT,'r')
    plt.title("{} on {:02d}-{:02d}-{:04d}".format(name,day,month,year))
    plt.xlabel('Azimuth [°]')
    plt.ylabel('Elevation [°]')
    plt.text(270,36,'Geo satellites',color='green') 
    plt.text(270,25,'Path sun',color='red') 
    plt.text(270,15,'Sagittarius',color='blue') 
    plt.text(270,43,'Telescope FOV',color='red') 
    plt.xticks(np.arange(0, 360+1, 30.0))
    plt.xlim(60,330)
    plt.ylim(0,50)
    plt.grid()
    for h in range(0,len(hh),3):
        if ((AZI[h]> AziMin) and  (AZI[h]<AziMax) and (ELE[h]>EleMin) and (ELE[h]<EleMax)):
            plt.text(AZI[h],ELE[h],hh[h])
            
    sun = ephem.Sun()
    for h in range(0,24):
        for m in range(0,60,10):
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year,month,day,h,m,0)
            sun.compute(MyHome)
            Sazi = math.degrees(sun.az)
            Sele = math.degrees(sun.alt)
            # if ((m%10)==0 and (h>6) and (h<16)):
            #     print("{:5.2f}, {:4.1f}, {:4.1f}".format(h+m/60,Sazi,Sele))
            plt.plot(Sazi,Sele,'.r')
        if ((h%1)==0):
            if ((Sazi> AziMin) and  (Sazi<AziMax) and (Sele>EleMin) and (Sele<EleMax)):
                plt.text(Sazi,Sele,h)
   
    for sat in range (0,360,5):
        a,e = GeoSat(sat)
        if (e > 0):
                plt.plot(a,e,'.g') 
        plt.plot(360-a,e,'.g') 

    plotfilename = '{:4d}-{:02d}-{:02d}_'.format(year,month,day) + name
    plt.savefig(plotfilename)

    plt.show()

#------------------------------------------------------------------------------
    
# MyAstro(83.5 , 22.0,'Taurus A A in Locarno Monti'     )    
# MyAstro(350.8, 58.8,'Cassiopeia A A in Locarno Monti' )    
# MyAstro(299.8, 40.7,'Cygnus A in Locarno Monti'     )    
# MyAstro(187.5, 12.4,'Virgo A (M87) A in Locarno Monti')    
# MyAstro(83.8 , -5.4,'Orion A A in Locarno Monti'      )    
MyAstro(265.5,-28.8,'Meteo Swiss in Locarno Monti') # Sagittarius plot 

#------------------------------------------------------------------------------
