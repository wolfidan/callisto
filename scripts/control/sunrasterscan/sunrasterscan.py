# -*- coding: utf-8 -*-
"""
This script does a raster scan of the sun, similarly to what the Leonardo DX50 weather radar is doing
Two operating modes are supported:
- coarse: scans the whole sky to find for the possible aziimuth and elevation of the sun
- fine: scans around the estimated aziimuth and elevation of the sun, to fine-tune the estimation of the coarse scan
"""
#%%
import datetime
import numpy as np
import time
import configparser
import time
import os
import atexit
from typing import Optional, Tuple

from callisto_lib.classes import SunEphem, PythonLogger, AntennaTracker
from callisto_lib.constants import DIR_CONFIG
from callisto_lib.math_utils import generate_scanning_pattern

# def make_sun_raster_scan(mode: str = "coarse", sun_pos: tuple[float, float] | None = None) -> None:
#------------------------------------------------------------------------------

mode = "coarse"
sun_pos = None

#%%
# Create python logger
logger = PythonLogger("sunrasterscan")

# Read configuration file
config = configparser.ConfigParser()
config.read(os.path.join(DIR_CONFIG, "config_sunrasterscan.ini"))

#------------------------------------------------------------------------------
# Tracker parameter
max_range  = config.getfloat('Tracker', 'max_range')
azi_min    = config.getfloat('Tracker', 'azi_min')
azi_max    = config.getfloat('Tracker', 'azi_max')
ele_min    = config.getfloat('Tracker', 'ele_min')
ele_max    = config.getfloat('Tracker', 'ele_max')
azi_park   = config.getfloat('Tracker', 'azi_park')
ele_park   = config.getfloat('Tracker', 'ele_park')
azi_ref    = config.getfloat('Tracker', 'azi_ref')
ele_ref    = config.getfloat('Tracker', 'ele_ref')
azi_dir    = config.getfloat('Tracker', 'azi_dir')
ele_dir    = config.getfloat('Tracker', 'ele_dir')
plane_corr = config.getfloat('Tracker', 'plane_corr')
com_port   = config.get('Tracker', 'com_port')

#%%
#------------------------------------------------------------------------------
# Create antenna tracker
now = datetime.datetime.now(datetime.timezone.utc) 
nowstr = now.strftime("%Y%m%dT%H%M")
diseqlog_file = f"DISEQ-{nowstr}-Sunraster-{mode}.txt"
tracker = AntennaTracker(com_port, 
                        azi_dir, azi_ref,
                        ele_dir, ele_ref,
                        plane_corr,
                        python_logger = logger,
                        diseqlog_file = diseqlog_file)

#------------------------------------------------------------------------------
# Scanning parameters
mode_cap = mode.title()
dwell_time = config.getfloat(mode_cap,'dwell_time')
transition_time = config.getfloat(mode_cap,'transition_time')
res_azi = config.getfloat(mode_cap,'res_azi')
res_ele = config.getfloat(mode_cap,'res_ele')
width_azi = config.getfloat(mode_cap,'width_azi')
width_ele =  config.getfloat(mode_cap,'width_ele')

if sun_pos:
    sun_azi = sun_pos[0]
    sun_ele = sun_pos[1]
else:
    # Create Sun Ephem
    sun_ephem = SunEphem(config.get('Location','Longitude'),
        config.get('Location','Latitude'),
        config.getfloat('Location','Elevation'),
        config.getfloat('Location','Temperature'),
        config.getfloat('Location','Pressure'))
    # Get it from ephemeridis
    sun_azi, sun_ele = sun_ephem.get_sun_angles(now)

# Get extent of scan    
extent_azi = (max(sun_azi - width_azi/2., azi_min),
                min(sun_azi + width_azi/2, azi_max))
extent_ele = (max(sun_ele - width_ele/2., ele_min),
                min(sun_ele + width_ele/2, ele_max))

atexit.register(tracker.close)  # Ensure port is closed on exit

# Create grid of azi,el
all_ele = np.arange(extent_ele[0], extent_ele[1]+res_ele, res_ele)
all_azi = np.arange(extent_azi[0], extent_azi[1]+res_azi, res_azi)

scanning_pattern = generate_scanning_pattern(all_azi, all_ele)
npoints = len(scanning_pattern)

logger.info(f"========================")
logger.info(f"Starting {mode} sun scan")
logger.info(f"Azi sun: {sun_azi:6.3f} deg")
logger.info(f"Ele sun: {sun_ele:6.3f} deg")
logger.info(f"Azi scan range: {extent_azi[0]:6.3f}-{extent_azi[1]:6.3f} deg, step = {res_azi} deg")
logger.info(f"Ele scan range: {extent_ele[0]:6.3f}-{extent_ele[1]:6.3f} deg, step = {res_ele} deg")
logger.info(f"Nb sampling points: {npoints}")
logger.info(f"Estimated scan time:  {(npoints*(dwell_time + transition_time))/60:4.3f} min")

npoints_by_percent = int(npoints // 100)
for i, (azi, ele) in enumerate(scanning_pattern):
    tracker.send_position(azi, ele, transition_time=transition_time/2.)
    if (i % npoints_by_percent) == 0:
        logger.info(f"{int((i / npoints)*100)}% done")
    time.sleep(dwell_time)

logger.info(f"Finished {mode} sun scan")
logger.info(f"Closing tracker")
tracker.close()
logger.info(f"========================")