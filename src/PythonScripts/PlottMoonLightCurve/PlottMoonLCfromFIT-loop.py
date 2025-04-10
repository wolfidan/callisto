# -*- coding: utf-8 -*-
"""
Created on Tuesday January 14, 2020
Script allows to read and plott several FIT-files from local folder.
Edit graphic-file extension if required png, pdf, jpg, ....

@author: Christian Monstein 2020-06-20,2023-07-27
"""
#-------------------------------------------------------------------------------

from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np
import glob

#-------------------------------------------------------------------------------

MyPath = '' # Edit in case the file is elsewhere, outside local folder, use '/' instead of '\'
MyFile = 'NASA_20230904_*_02.fit'
paths = glob.glob(MyPath + MyFile) # search the folder
paths.sort() # force alphabetical order of the files
print ('\nFIT-files found: ',(len(paths)),paths)

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
    dT         = hdu[0].header['CDELT1']     # extract time resolution
    datum      = hdu[0].header['DATE-OBS']   # take first file
    T0         = hdu[0].header['TIME-OBS']   # take first file
    instrument = hdu[0].header['INSTRUME']
    content    = hdu[0].header['CONTENT']
    # https://docs.astropy.org/en/stable/io/fits/
    # print (hdu.info())
    # print (hdu[0].header)
    # print (hdu[1].header)
    
#---------------------------------------------------------------------------------

for path in paths:
    data   = get_data_from_fits(path)
    Tstart = get_T0_from_fits(path)
    st = Tstart[0:2] + Tstart[3:5]
    
    dB = data/255*2500.0/25.4 # raw-date -> dB
    Dlin = 10.0**(dB/10) # dB -> linear power
    cols = len(Dlin[0])
    T = np.arange(0,cols,1) * dT # sampling rate usually 2 Hz
    
    #LC = Dlin[92,:] # average of all but only 8 bit
    LC = np.mean(Dlin[40:80,:], axis=0) # take only frequencies without SAT-TV / 0.25 sec FIT-files
    
    plt.figure(figsize=(16,5))# present each 15 minute file
    plt.step(T,LC,'-',color='green',label=st)
    plt.xlabel("Relative time [sec]")
    plt.ylabel("Intensity [digits]")
    plt.yscale('log')
    plt.title('Light curve raw data of: '+path)
    plt.legend(loc='lower right')
    plt.savefig("XbandRasterFlyScan_LC_%s_%s.png"%(datum.replace("/", "_"),path),bbox_inches="tight")

#-------------------------------------------------------------------------------
