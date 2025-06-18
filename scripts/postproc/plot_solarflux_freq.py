# -*- coding: utf-8 -*-

description = """
Description:

    Plot solarflux measurements generated from FIT-files. This is an extension of 
    Philipp's script "plot_solarflux.py." This version can handle also the 
    spectrogram plots at all frequencies.

    @author: Andrea F. Battaglia

    History:
        - 2024/12/19 [Andrea F. Battaglia]: first version

"""

import os
import datetime
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import astropy.units as u
import pickle
import numpy as np

from astropy.time import Time
from tools.utils import checkdirs
from sunpy.net import Fido, attrs as a
from datetime import timedelta
from sunpy import timeseries as ts
from matplotlib.dates import DateFormatter
from sunpy.time import parse_time
from matplotlib.colors import LogNorm

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import warnings
warnings.filterwarnings('ignore')

def get_datetime(date, time):
    dt = datetime.datetime(int(date[0:4]), int(date[5:7]), int(date[8:10]), int(time[0:2]), int(time[3:5]))
    return dt

def main():

    # Adding argparser allows to add arguments when running the code from the
    # terminal. Below, two arguments are added:
    #   - for the selection of the day
    #   - to store the plot somewhere
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=description, epilog='')
    
    parser.add_argument(
        '-day', type=str, metavar=str.__name__, 
        action='store', required=False, 
        help='Choose day to analyze, e.g. 20240510 (default: today)')
    
    parser.add_argument(
        '-store', action='store_true', required=False, default=False,
        help='whether to store figures or not (default: %(default))')
    
    parser.add_argument(
        '-remove_pkl', action='store_true', required=False, default=False,
        help='whether to remove the original pkl file (default: False)')
    
    args = parser.parse_args()
    day = args.day
    store = args.store
    remove_pkl = args.remove_pkl

    ### Select three frequencies for plotting the integrated time profiles
    ## in MHz
    #df = 7.62304688 # nominal frequency resolution
    fr1 = 11065
    fr2 = 10800
    fr3 = 10640

    directory_procdata = 'C:\\xrt\\output\\data\\proc'
    directory_fig = 'C:\\xrt\\output\\fig\\solarflux'
    checkdirs([directory_procdata, directory_fig])

    ## Find the pkl file
    if not day:
        day_datetime = datetime.date.today()
        day = day_datetime.strftime("%Y%m%d")
    filepath_pkl = os.path.join(directory_procdata, day, 'processed-obs_')+day+'.pkl'
    if not os.path.exists(filepath_pkl):
        print('No such file: %s' % filepath_pkl)
    
    ## Define the path of the figure to store
    filepath_outputfig = os.path.join(directory_fig, "callisto_meteoswiss_%s.png" % day)

    
    with open(filepath_pkl, 'rb') as pkl:
        date = pickle.load(pkl)
        raw_data = pickle.load(pkl)
        reduced_data = pickle.load(pkl)
        freq_axis_MHz = pickle.load(pkl)
        time_axis_raw = pickle.load(pkl)
        solar_obs = pickle.load(pkl)
        bkg_obs = pickle.load(pkl)
        Thot_linpow = pickle.load(pkl)
        Tcold_linpow = pickle.load(pkl)
        start_times = pickle.load(pkl)
        time_axis_sun = pickle.load(pkl)
        time_axis_Thot = pickle.load(pkl)
        time_axis_Tcold = pickle.load(pkl)
        Tcold = pickle.load(pkl)
        Thot = pickle.load(pkl)

    ntimes = len(time_axis_raw)
    
    ## here we download the solar X-ray data, to check the solar activity
    # define where to store the downloaded data
    goes_data_path = 'C:\\xrt\\output\\data\\goes'

    now_utc = datetime.datetime.utcnow()
    if now_utc - pd.DateOffset(days=6) > time_axis_raw[ntimes-1]:
        
        # download the GOES/XRS data via the standard pipeline
        result_xrs = Fido.search(a.Time(str(time_axis_raw[0]), str(Time(time_axis_raw[ntimes-1]) + timedelta(hours=2))), 
                                a.Instrument("XRS"), a.Resolution("avg1m"))
        file_goes = Fido.fetch(result_xrs[0, 0:2], path=goes_data_path)
        
        # create the timeseries structure and truncate it to the time of interest
        goes_xrs = ts.TimeSeries(file_goes, concatenate=True)
        
    else:
        
        # Since the requested data are too close to the observing day, we
        # need to request nrt (near-real-time data), which requires a different coding
        goes_json_data = pd.read_json("https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json")
        goes_short = goes_json_data[goes_json_data["energy"] == "0.05-0.4nm"]
        goes_long = goes_json_data[goes_json_data["energy"] == "0.1-0.8nm"]
        time_array = parse_time(goes_short["time_tag"])
        units = dict([("xrsa", u.W/u.m**2), ("xrsb", u.W/u.m**2)])
        meta = dict({"instrument": "GOES X-ray sensor", "measurements": "primary", "type": "quicklook"})
        goes_data = pd.DataFrame({"xrsa": goes_short["flux"].values, "xrsb": goes_long["flux"].values}, index=time_array.datetime)
        goes_xrs = ts.TimeSeries(goes_data, meta, units, source="xrs")

    goes_xrs = goes_xrs.truncate(time_axis_raw[0], time_axis_raw[ntimes-1])

    ## Calculate the time profiles by averaging the fluxes
    idx_f1 = np.abs(freq_axis_MHz - fr1).argmin()
    idx_f2 = np.abs(freq_axis_MHz - fr2).argmin()
    idx_f3 = np.abs(freq_axis_MHz - fr3).argmin()
    lc_sun_f1 = solar_obs[idx_f1,:]
    lc_sun_f2 = solar_obs[idx_f2,:]
    lc_sun_f3 = solar_obs[idx_f3,:]
    lc_Thot_f1 = Thot_linpow[idx_f1,:]
    lc_Thot_f2 = Thot_linpow[idx_f2,:]
    lc_Thot_f3 = Thot_linpow[idx_f3,:]
    lc_Tcold_f1 = Tcold_linpow[idx_f1,:]
    lc_Tcold_f2 = Tcold_linpow[idx_f2,:]
    lc_Tcold_f3 = Tcold_linpow[idx_f3,:]
    
    Yfact1 = lc_Thot_f1 / lc_Tcold_f1
    Yfact1 = 10.0*np.log10(Yfact1)
    Yfact2 = lc_Thot_f2 / lc_Tcold_f2
    Yfact2 = 10.0*np.log10(Yfact2)
    Yfact3 = lc_Thot_f3 / lc_Tcold_f3
    Yfact3 = 10.0*np.log10(Yfact3)

    dt_axis = pd.to_timedelta(2, unit='minutes')

    #****************************** PLOT ******************************

    fig = plt.figure(figsize=(12,6.5))
    plt.rcParams.update({'font.size': 8})

    # GOES/XRS X-ray subplot
    axg = fig.add_subplot(321)
    goes_xrs.plot(axes=axg)
    axg.set_title('Solar activity - GOES/XRS X-ray flux for '+ time_axis_raw[0].strftime('%Y-%m-%d'))
    axg.xaxis_date()  # Use date format for the x-axis
    axg.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display

    # Plot the raw data
    axr = fig.add_subplot(322)
    #extent = timeax[0], timeax[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
    #im = ax1.imshow(tmp_solar_obs, aspect='auto', vmin=100, vmax=600, extent=extent)
    extent = time_axis_raw[0], time_axis_raw[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
    imr = axr.imshow(raw_data, aspect='auto', origin='lower', extent=extent, norm=LogNorm(vmin=1e4, vmax=0.5e7))
    caxr = fig.add_axes([axr.get_position().x0 + 0.05,
                        axr.get_position().y0 + 0.305,
                        axr.get_position().width,
                            0.02])
    cbar = fig.colorbar(imr, cax=caxr, orientation='horizontal')
    label = cbar.set_label('Raw data; linear power', color='black')
    label_text = cbar.ax.get_xaxis().get_label()
    label_text.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))
    cbar.ax.tick_params(axis='x', which='both', direction='in', pad=-8.5, labelsize=6)
    for label in cbar.ax.xaxis.get_ticklabels():
        label.set_verticalalignment('top')
    axr.set_xlabel('UTC')
    axr.set_ylabel('Freq [GHz]')
    #axr.set_title('Raw data')
    axr.xaxis_date()  # Use date format for the x-axis
    axr.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display
    axr.hlines(fr1/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)
    axr.hlines(fr2/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)
    axr.hlines(fr3/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)

    # Plot the solar time profiles
    axlc = fig.add_subplot(323)
    axlc.plot(time_axis_sun, lc_sun_f1, '.', alpha=0.8, label=str(fr1/1e3)+' GHz')
    axlc.plot(time_axis_sun, lc_sun_f2, '.', alpha=0.8, label=str(fr2/1e3)+' GHz')
    axlc.plot(time_axis_sun, lc_sun_f3, '.', alpha=0.8, label=str(fr3/1e3)+' GHz')
    axlc.set_ylabel('CALLISTO flux density [sfu]')
    axlc.set_xlabel('UTC')
    #axlc.legend()
    axlc.legend(bbox_to_anchor=(0., 1.02, 1, 0.2), loc='lower left', 
          mode="expand", borderaxespad=0, ncol=3)
    axlc.set_yscale('log')
    #axlc.set_title('Callisto@MeteoSwiss solar flux for '+ time_axis_raw[0].strftime('%Y-%m-%d'))
    axlc.xaxis_date()  # Use date format for the x-axis
    axlc.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display

    # Plot the reduced observations in sfu
    axred = fig.add_subplot(324)
    extent = time_axis_raw[0], time_axis_raw[-1], freq_axis_MHz[0]/1e3, freq_axis_MHz[-1]/1e3
    imred = axred.imshow(reduced_data, aspect='auto', origin='lower', extent=extent, vmin=240, vmax=610)
    caxrred = fig.add_axes([axred.get_position().x0 + 0.05,
                            axred.get_position().y0 + 0.255,
                            axred.get_position().width,
                            0.02])
    cbarred = fig.colorbar(imred, cax=caxrred, orientation='horizontal')
    cbarred.set_label('Reduced data; Solar flux [sfu]', color='white')
    label_text = cbarred.ax.get_xaxis().get_label()
    label_text.set_bbox(dict(facecolor='black', alpha=0.5, edgecolor='none'))
    cbarred.ax.tick_params(axis='x', direction='in', pad=-8.5, labelsize=6)
    for label in cbarred.ax.xaxis.get_ticklabels():
        label.set_verticalalignment('top')
    axred.set_xlabel('UTC')
    axred.set_ylabel('Freq [GHz]')
    #axred.set_title('Reduced data')
    axred.xaxis_date()  # Use date format for the x-axis
    axred.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display
    axred.hlines(fr1/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)
    axred.hlines(fr2/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)
    axred.hlines(fr3/1e3, time_axis_raw[0], time_axis_raw[-1], colors='red', alpha=0.5)
    
    # Plot the Thot and Tcold
    axT = fig.add_subplot(325)
    axT.plot(time_axis_Thot+dt_axis, lc_Thot_f1, '.', c='brown', alpha=0.8, label='Thot, '+str(fr1/1e3)+' GHz')
    axT.plot(time_axis_Thot, lc_Thot_f2, '.', c='red', alpha=0.8, label='Thot, '+str(fr2/1e3)+' GHz')
    axT.plot(time_axis_Thot-dt_axis, lc_Thot_f3, '.', c='orange', alpha=0.8, label='Thot, '+str(fr3/1e3)+' GHz')
    axT.plot(time_axis_Tcold+dt_axis, lc_Tcold_f1, '.', c='green', alpha=0.8, label='Tcold, '+str(fr1/1e3)+' GHz')
    axT.plot(time_axis_Tcold, lc_Tcold_f2, '.', c='blue', alpha=0.8, label='Tcold, '+str(fr2/1e3)+' GHz')
    axT.plot(time_axis_Tcold-dt_axis, lc_Tcold_f3, '.', c='cyan', alpha=0.8, label='Tcold, '+str(fr3/1e3)+' GHz')
    axT.set_ylabel('Linear power')
    axT.set_xlabel('UTC')
    #axT.legend()
    axT.legend(bbox_to_anchor=(0.05, 1.02, 0.95, 0.2), loc='lower left', 
          mode="expand", borderaxespad=0, ncol=3)
    axT.xaxis_date()  # Use date format for the x-axis
    axT.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display

    # Plot the Y factor
    axY = fig.add_subplot(326)
    axY.plot(time_axis_Thot+dt_axis, Yfact1, '.', c='brown', alpha=0.8, label=str(fr1/1e3)+' GHz')
    axY.plot(time_axis_Thot, Yfact2, '.', c='gray', alpha=0.8, label=str(fr2/1e3)+' GHz')
    axY.plot(time_axis_Thot-dt_axis, Yfact3, '.', c='black', alpha=0.8, label=str(fr3/1e3)+' GHz')
    axY.set_ylabel('Y-factor [dB]')
    axY.set_xlabel('UTC')
    #axY.legend()
    axY.legend(bbox_to_anchor=(0.00, 1.02, 1, 0.2), loc='lower left', 
          mode="expand", borderaxespad=0, ncol=3)
    axY.xaxis_date()  # Use date format for the x-axis
    axY.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display

    plt.tight_layout()

    if remove_pkl:
        os.remove(filepath_pkl)
        store = True # if we remove the file, by default we store the image, else it does not make snse to call this function

    if store:
        print('Export figure: %s' % filepath_outputfig)
        plt.savefig(filepath_outputfig,bbox_inches="tight")    
    
    if not store:
        plt.show()
    #plt.title('X-band solar radio flux at {:7.4f} GHz with {:6.1f} SFU'.format(f/1e9,Smed))    
    #plt.savefig("XbandRasterFlyScan_LC_%s.png"%(datum.replace("/", "_")),bbox_inches="tight")

    return

if __name__ == '__main__':
    main()