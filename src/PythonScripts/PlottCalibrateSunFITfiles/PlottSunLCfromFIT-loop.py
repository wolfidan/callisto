# -*- coding: utf-8 -*-
"""
Created on Tuesday January 14, 2020
Script allows to read and plott several FIT-files from local folder.
Edit graphic-file extension if required png, pdf, jpg, ....

@author: Christian Monstein 2020-06-20,2023-07-27,2023-09-12
"""
#-------------------------------------------------------------------------------

from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np
import scipy.constants as CON
import glob
import sys
import os

#-------------------------------------------------------------------------------

debug = True # True=all FIT-files -> plot, False=only result
Tcold = 12 # from literature
Thot = 273.15 + 12.0 # from outdoor sensor
Terr = 0 # timing error FIT-files in seconds, should be 0 in normal cases
LearmonthFlux = 343.8 # check here: https://owenduffy.net/calc/qsrf/index.htm
YcalMin = 3.0 # minimum Y-factor during calibration to accept result; quality measure

# C:\\xrt\\src\\PythonScripts\\PlottCalibrateSunFITfiles\\
#MyPath = 'C:\\xrt\\src\\PythonScripts\\PlottCalibrateSunFITfiles\\' # Edit in case the file is elsewhere, outside local folder, use '/' instead of '\'
#MyPath = ''
MyPath = 'C:/xrt/output/data/raw/FITfiles'
MyFile = 'meteoswiss_20240510_*_01.fit'
#MyFile = 'SWISS-METEO_20231127_*_01.fit'
print(MyPath)
#print(MyPath+MyFile)
paths = glob.glob(os.path.join(MyPath,MyFile)) # search the folder
paths.sort() # force alphabetical order of the files
print('\nFIT-files found: ',(len(paths)),paths)

#-------------------------------------------------------------------------------

def get_data_from_fits(path):
    with fits.open(path) as hdu:
        data = hdu[0].data.astype(np.uint8) # float16, float32, uin8 depends on PC
    return data

#-------------------------------------------------------------------------------

def get_T0_from_fits(path):
    with fits.open(path) as hdu:
        T0 = hdu[0].header['TIME-OBS'] 
    return T0

#-------------------------------------------------------------------------------

with fits.open(paths[0]) as hdu:
    data       = hdu[0].data
    freqs      = hdu[1].data  ['Frequency'][0] # extract frequency axis
    timeax     = hdu[1].data  ['Time'][0]      # extract time axis
    dT         = hdu[0].header['CDELT1']       # extract time resolution
    datum      = hdu[0].header['DATE-OBS']     # take first file
    T0         = hdu[0].header['TIME-OBS']     # take first file
    instrument = hdu[0].header['INSTRUME']
    content    = hdu[0].header['CONTENT']
    # https://docs.astropy.org/en/stable/io/fits/
    # print (hdu.info())
    # print (hdu[0].header)
    # print (hdu[1].header)

   
# print(paths)
# print(hdu.info())
# print(hdu[0].header)
# print(hdu[1].header)
# print(data)
# sys.exit()
    
#---------------------------------------------------------------------------------

Flux = [] # Solar flux list
X    = [] # Decimal time of FIT-file

print(paths)
for path in paths:
    data   = get_data_from_fits(path)
    Tstart = get_T0_from_fits(path)
    st = Tstart[0:2] + Tstart[3:5]
    x = int(Tstart[0:2]) + int(Tstart[3:5])/60 # Decimal hours
    #print ('Data shape: ',data.shape)
    
    #------------------------------------------------------------------------------
    dB = data/255*2500.0/25.4 # raw-date -> dB
    Dlin = 10.0**(dB/10) # dB -> linear power
    cols = len(Dlin[0])
    T = np.arange(0,cols,1) * dT # sampling rate usually 2 Hz
    
    # LC = np.mean(Dlin[140:180,:], axis=0) # take only frequencies without SAT-TV / 0.5 sec FIT-files with 200 channels
    #LC = np.mean(Dlin[40:80,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files with 100 channels
    # Plot spectrum to find pointers, here 40:80
    LC = np.mean(Dlin[81:86,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files with 100 channels
    
    #-------------------------------------------------------------------------------
    
    print ('\n\nResult 15 minute observation: ',datum,' at ',Tstart)
    T1 = int((10-Terr)/dT) # start Icold, sec -> pixel
    T2 = int((60-Terr)/dT)
    Icold = np.mean(LC[T1:T2])
    
    T3 = int(( 80-Terr)/dT) # start Ihot, sec -> pixel
    T4 = int((120-Terr)/dT)
    Ihot  = np.mean(LC[T3:T4])
    
    Ycal  = Ihot/Icold
    YcaldB = 10.0*np.log10(Ycal)
    print ('Icold {:8.1f} digit, Ihot {:8.1f} digit, YcaldB {:5.2f} dB '.format(Icold,Ihot,YcaldB))
    
    s = np.std(LC[T1:T2] )
    snr = (Ihot-Icold)/s
    #print ('snr {:6.1f} = {:5.1f} dB'.format(snr,10.0*np.log10(snr)))
    
    T5 = int((170-Terr)/dT) # start Isunscan, sec -> pixel
    T6 = int((900-Terr)/dT) # 900
    sun = LC[T5:T6]
    sunpeak = np.max(sun)
    #print ('Peak sun {:8.1f} digit '.format(sunpeak))
    
    i = 0
    I = 0
    MyFluxLimit = 0.97 # take only values above 97% of peak flux
    for isun in sun:
        if (isun > MyFluxLimit*sunpeak):
            I = I + isun
            i = i + 1
    Isun = I/i
    #print ('Average sun (0.97*SunPeak) {:8.1f} digit from {:4.0f} measured pixels'.format(Isun,i))
    
    Tsun = (Thot - Tcold) / (Ihot - Icold) * (Isun - Icold) + Tcold
    #print ('\nAntenna temperature sun: {:6.1f} Kelvin'.format(Tsun))
    Ysun = (Isun/Icold)
    YsundB =  10.0*np.log10(Ysun)
    #print ('Ysun {:5.2f} dB '.format(YsundB))
    
    gain = 36 # # according to data sheet of TRIAX TDS 65 A
    G    = 10.0**(gain/10.) 
    sfu  = 1e22
    f    = 11059e6 # from line ~76 in case of 0.25 sec FIT-file
    lamb = CON.c/f
    Aeff = lamb**2 * G / (4*np.pi)
    #print ('Effecticve aperture: {:5.3f} m^2'.format(Aeff))
    
    Ssun = 2*CON.Boltzmann*Tsun/Aeff*sfu
    print ('Solar radio flux: {:6.1f} SFU'.format(Ssun))

    Sback = (2*CON.Boltzmann*Tcold)/(lamb**2*G)*4*np.pi * sfu
    print ('Background flux:            {:6.1f} SFU'.format(Sback))
    print ('Final flux above background:{:6.1f} SFU'.format(Ssun-Sback))
    
    if (YcaldB > YcalMin): # Check quality of calibration Thot/Tcold > limit
        X.append(x)  # Start-time FIT-file
        Flux.append(Ssun)
        
        if (debug):
            st = 'Flux {:6.1f} SFU at  {:6.3f} GHz'.format(Ssun,f/1e9)
            plt.figure(figsize=(16,5))# present each 15 minute file
            plt.step(T,LC,'-',color='green',label=st)
            plt.xlabel("Relative time [sec]")
            plt.ylabel("Intensity [digits]")
            plt.yscale('log')
            plt.title('Light curve raw data of: '+path)
            plt.legend(loc='lower right')
            plt.savefig("XbandRasterFlyScan_LC_%s_%s.png"%(datum.replace("/", "_"),path),bbox_inches="tight")


#-------------------------------------------------------------------------------

if len(Flux)>0:
    Smed = np.median(Flux)
    Srms = np.std(Flux)
else:
    Smed = np.nan
    Srms = np.nan

print ('\nFinal results:')
avg = '{:6.1f} SFU'.format(Smed)
print ('Solar radio flux median:    {:6.1f} SFU'.format(Smed))
print ('Solar flux OwenDuffy.net:   {:6.1f} SFU'.format(356.4))  # reference flux
print ('Solar radio flux std:       {:6.1f} SFU'.format(Srms))

# Background flux estimation to be subtracted

Sback = (2*CON.Boltzmann*Tcold)/(lamb**2*G)*4*np.pi * sfu
print ('Background flux:            {:6.1f} SFU'.format(Sback))
print ('Final flux above background:{:6.1f} SFU (median)'.format(Smed-Sback))

plt.figure(figsize=(8,5)) # Flux versus time
plt.bar(X,Flux,color='green',width=0.2)
plt.xlabel("Decimal time [UT] on " + datum)
plt.ylabel("Flux [sfu]")
plt.grid('both')
plt.title('X-band solar radio flux at {:7.4f} GHz with {:6.1f} SFU'.format(f/1e9,Smed))
plt.savefig("XbandRasterFlyScan_LC_%s.png"%(datum.replace("/", "_")),bbox_inches="tight")

#-------------------------------------------------------------------------------
