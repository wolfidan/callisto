# -*- coding: utf-8 -*-

description = """
Description:

    Plot solarflux measurements generated from FIT-files

    @author: Philipp Schmid

    2024/05/10: make draft based on Christian Monstein's script "PlottSunLCfromFIT-loop.py"

    #### compute Aeff at each frequency
    #### define frequencies over which to integrate
    #### ...
    
    History:
        - 2024/12/11 [Andrea F. Battaglia]: added the argument -goes, to plot the
                                            solar X-ray flux (to check the solar activity).
                                            To activate this feature, add '-goes 1'
        - 2024/12/13 [Andrea F. Battaglia]: added the argument -raw, in order to plot the 
                                            raw data. This is helpful to check any systematic
                                            uncertainty. To activate this feature, add '-raw 1'

"""

import os
import datetime
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import astropy.units as u
import pickle

from astropy.time import Time
from tools.utils import checkdirs
from sunpy.net import Fido, attrs as a
from datetime import timedelta
from sunpy import timeseries as ts
from matplotlib.dates import DateFormatter
from sunpy.time import parse_time

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
        action='store', required=True, 
        help='Choose day to analyze, e.g. 20240510')
    
    parser.add_argument(
        '-store', action='store_true', required=False, default=False,
        help='whether to store figures or not (default: %(default))')
    
    parser.add_argument(
        '-goes', required=False, default=False, 
        help='Choose whether to plot the solar X-ray flux (GOES/XRS)'
    )

    parser.add_argument(
        '-raw', required=False, default=False, 
        help='Choose whether to plot the raw linear power curves'
    )

    args = parser.parse_args()
    day = args.day
    store = args.store
    goes = args.goes
    raw = args.raw

    directory_procdata = 'C:\\xrt\\output\\data\\proc'
    directory_fig = 'C:\\xrt\\output\\fig\\solarflux'
    checkdirs([directory_procdata, directory_fig])

    ## Find the files
    filepath_procdata = os.path.join(directory_procdata, day, 'solarflux.csv')
    if not os.path.exists(filepath_procdata):
        print('No such file: %s' % filepath_procdata)
    filepath_outputfig = os.path.join(directory_fig, "solarflux_%s.pdf" % day)

    ## If the argument "-raw" is set, then open the pkl file
    if raw:
        filepath_pkl = os.path.join(directory_procdata, day, 'time-profiles_sun_Thot_Tcold_')+day+'.pkl'
        if not os.path.exists(filepath_pkl):
            print('No such file: %s' % filepath_pkl)

        with open(filepath_pkl, 'rb') as pkl:
            date = pickle.load(pkl)
            raw_data = pickle.load(pkl)
            freq_axis_MHz = pickle.load(pkl)
            time_axis_raw_s = pickle.load(pkl)
            start_times = pickle.load(pkl)
            lc_sun_linpow = pickle.load(pkl)
            lc_Thot_linpow = pickle.load(pkl)
            lc_Tcold_linpow = pickle.load(pkl)
            time_axis_sun_s = pickle.load(pkl)
            time_axis_Thot_s = pickle.load(pkl)
            time_axis_Tcold_s = pickle.load(pkl)


    with open(os.path.join(filepath_procdata), 'r') as f:
        Tcold = f.readline().strip().split(': ')[1]
        Thot = f.readline().strip().split(': ')[1]
        freq_MHz = f.readline().strip().split(': ')[1]
        dfreq_MHz = f.readline().strip().split(': ')[1]

    df = pd.read_csv(filepath_procdata,delimiter=', ', engine='python', comment='#', skiprows=4)
    df['datetime'] = df.apply(lambda dp: get_datetime(dp['date'], dp['time_start']), axis=1)
    datetime_callisto = pd.to_datetime(df['datetime'])
    ntimes = len(datetime_callisto)
    
    # here we download the solar X-ray data, to check the solar activity
    if goes == '1':

        # define where to store the downloaded data
        goes_data_path = 'C:\\xrt\\output\\data\\goes'

        now_utc = datetime.datetime.utcnow()
        if now_utc - pd.DateOffset(days=6) > datetime_callisto[ntimes-1]:
            
            # download the GOES/XRS data via the standard pipeline
            result_xrs = Fido.search(a.Time(str(datetime_callisto[0]), str(Time(datetime_callisto[ntimes-1]) + timedelta(hours=2))), 
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

        goes_xrs = goes_xrs.truncate(datetime_callisto[0], datetime_callisto[ntimes-1])

    min_y = min(df['Ssun_minus_Sback_sfu']) / 2
    if min_y <= 0: min_y = 100
    # distinguish the two plots. If the -goes argument is called, then do a different plot
    if goes == '1':

        fig = plt.figure(figsize=(9,6))

        # GOES/XRS X-ray subplot
        axg = fig.add_subplot(211)
        goes_xrs.plot(axes=axg)
        axg.set_title('Solar activity - GOES/XRS X-ray flux for '+ datetime_callisto[0].strftime('%Y-%m-%d'))
        axg.xaxis_date()  # Use date format for the x-axis
        axg.xaxis.set_major_formatter(DateFormatter("%H:%M"))  # Format time display

        # CALLISTO subplot
        time_axis_bar = pd.to_datetime(df['datetime']) + pd.to_timedelta(60*7.5, unit='s')
        ax = fig.add_subplot(212)#, sharex=axg)
        ax.bar(time_axis_bar, df['Ssun_minus_Sback_sfu'],color='green',width=datetime.timedelta(minutes=10))
        #ax.plot(df['datetime'], df['Ssun_minus_Sback_sfu'], '-o')
        ax.set_xlabel("start time of 15 min measurement interval")
        ax.set_ylabel("solarflux (noise subtracted) [sfu]")
        ax.grid('both')
        ax.set_ylim([min_y, max(df['Ssun_minus_Sback_sfu'])*2])
        ax.set_yscale('log')
        ax.set_title('Solar flux on day %s (Tcold: %s K, Thot: %s K, freq: %s MHz, dfreq: %s MHz)' % (
            day, Tcold, Thot, freq_MHz, dfreq_MHz), fontsize='small')

        myFmt = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)

    else:
            
        # fig, ax = plt.subplots()
        # ax.plot(df['datetime'], df['Ssun_minus_Sback_sfu'])
        # plt.show()

        time_axis_bar = pd.to_datetime(df['datetime']) + pd.to_timedelta(60*7.5, unit='s')
        fig, ax = plt.subplots(figsize=(7,4))
        ax.bar(time_axis_bar, df['Ssun_minus_Sback_sfu'],color='green',width=datetime.timedelta(minutes=10))
        #ax.plot(df['datetime'], df['Ssun_minus_Sback_sfu'], '-o')
        ax.set_xlabel("start time of 15 min measurement interval")
        ax.set_ylabel("solarflux [sfu]")
        ax.grid('both')
        ax.set_yscale('log')
        ax.set_ylim([min_y, max(df['Ssun_minus_Sback_sfu'])*2])
        ax.set_title('Solar flux on day %s (Tcold: %s K, Thot: %s K, freq: %s MHz, dfreq: %s MHz)' % (
            day, Tcold, Thot, freq_MHz, dfreq_MHz), fontsize='small')

        myFmt = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)

    
    ## If the argument raw is activated, plot the raw lightcurves
    if raw:
        #print(date)
        #print(start_times[0][0:8])
        this_start_time = pd.to_datetime(date+' '+start_times[0][0:8])
        time_axis_sun = this_start_time + pd.to_timedelta(time_axis_sun_s, unit='s')
        time_axis_Thot = this_start_time + pd.to_timedelta(time_axis_Thot_s, unit='s')
        time_axis_Tcold = this_start_time + pd.to_timedelta(time_axis_Tcold_s, unit='s')
        ax2 = ax.twinx()
        if raw == 'sun' or raw == 'all': ax2.plot(time_axis_sun, lc_sun_linpow, '.', c='tab:orange', markersize=2, label='sun')
        if raw == 'Thot' or raw == 'all': ax2.plot(time_axis_Thot, lc_Thot_linpow, '.', c='red', markersize=2, label='Thot')
        if raw == 'Tcold' or raw == 'all': ax2.plot(time_axis_Tcold, lc_Tcold_linpow, '.', c='blue', markersize=2, label='Tcold')
        ax2.set_ylabel('Linear power')
        ax2.legend()


    if store:
        print('Export figure: %s' % filepath_outputfig)
        plt.savefig(filepath_outputfig,bbox_inches="tight")

    plt.show()
    #plt.title('X-band solar radio flux at {:7.4f} GHz with {:6.1f} SFU'.format(f/1e9,Smed))    
    #plt.savefig("XbandRasterFlyScan_LC_%s.png"%(datum.replace("/", "_")),bbox_inches="tight")

    return

if __name__ == '__main__':
    main()