# -*- coding: utf-8 -*-
"""
Created on 2020-04-13, 2020-05-06, 2020-12-10
Azimuth and elevation of any fixstar, defined by right-ascension and declination at given date
@author: Christian Monstein
"""
#------------------------------------------------------------------------------

import ephem
import numpy as np
from matplotlib import pyplot as plt
import math
import os

#------------------------------------------------------------------------------
# Date to be chosen by the user
year  = 2023
month =  9
day   = 18
#------------------------------------------------------------------------------
# Location parameter for Bleien

MyHome = ephem.Observer()
MyHome.lon      = '-76.8423'  # east +°  GSFC
MyHome.lat      = '38.99149'  # north +°
MyHome.elev     = 48          # altitude in m asl
MyHome.temp     = 10          # °C
MyHome.pressure = 900         # mbar

#------------------------------------------------------------------------------
# Telescope parameter for dishes 5m and 7m

AziMin    =  180-77
AziMax    =  180+77
EleMin    =  -5 # 5m safety
EleMax    =  85.0
xT = [AziMin,AziMax,AziMax,AziMin,AziMin] # telescope frame azimuth
yT = [EleMin,EleMin,EleMax,EleMax,EleMin] # telescope frame elevation

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
        for m in range(0,60,10):
            #dt = datetime.datetime.now()
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year,month,day,h,m,0)
            
            fb.compute(MyHome);
            alt = np.rad2deg(1.*fb.alt)
            azi = np.rad2deg(1.*fb.az);
            hm = h+m/60
            HM.append(hm)
            AZI.append(azi)
            ELE.append(alt)
            st = '{:4.1f}'.format(hm)
            hh.append(st)

    
    plt.figure()
    plt.plot(AZI,ELE,'.b')
    plt.plot(xT,yT,'g')
    plt.title("{} on {:02d}-{:02d}-{:04d}".format(name,day,month,year))
    plt.xlabel('Azimuth [°]')
    plt.ylabel('Elevation [°]')
    plt.xticks(np.arange(0, 360+1, 30.0))
    plt.xlim(0,360)
    plt.ylim(0,90)
    for h in range(0,len(hh),3):
        if ((AZI[h]> AziMin) and  (AZI[h]<AziMax) and (ELE[h]>EleMin) and (ELE[h]<EleMax)):
            plt.text(AZI[h],ELE[h],hh[h])
            
    sun = ephem.Sun()
    for h in range(0,24):
        for m in range(0,60,10):
            #dt = datetime.datetime.now()
            MyHome.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
                    year,month,day,h,m,0)
            sun.compute(MyHome)
            Sazi = math.degrees(sun.az)
            Sele = math.degrees(sun.alt)
            if ((m%10)==0 and (h>6) and (h<16)):
                print("{:5.2f}, {:4.1f}, {:4.1f}".format(h+m/60,Sazi,Sele))
            plt.plot(Sazi,Sele,'.r')
            #print("Sun data: Azimuth   =%6.2f  Elevation   =%6.2f" % (azi-180.0,ele)) 
        if ((h%1)==0):
            if ((Sazi> AziMin) and  (Sazi<AziMax) and (Sele>EleMin) and (Sele<EleMax)):
                plt.text(Sazi,Sele,h)

    plotfilename = '{:4d}-{:02d}-{:02d}_'.format(year,month,day) + name
    plt.savefig(plotfilename)
    plt.show()

#------------------------------------------------------------------------------
try:
    os.remove("C:/Users/cmons/Documents/MyPython/AstroPlan/scheduler.txt")
except:
    print ("File not there\n")
    
# MyAstro(83.5 , 22.0,'Taurus A'     )    
#MyAstro(350.8, 58.8,'Cassiopeia A' )    
MyAstro(299.8, 40.7,'Cygnus A'     )    
# MyAstro(187.5, 12.4,'Virgo A (M87)')    
# MyAstro(83.8 , -5.4,'Orion A'      )    
# MyAstro(265.5,-28.8,'Sagittarius A')    

#------------------------------------------------------------------------------
