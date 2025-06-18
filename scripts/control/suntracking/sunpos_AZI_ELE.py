# -*- coding: utf-8 -*-
"""
This script tracks the Sun in azimuthal mode (azimuth and elevation)
Carefully check parameter in external file config.ini
It is sufficient to execute this script once every minute
Created on Thu Sep 14 20:13:50 2017
Updated: 09.05.2020, 05.05.2023, 13.05.2023, 04.08.2023, 23.08.2023
Updated; 23.12.2024 fixed positions for Tcold and Thot
Updated [Andrea]: 18.03.2025 added the scan of the sun for assessing the pointing offset
Updated [Daniel]: 20.05.2025 complete rewrite for script to run forever and with more options
@author: Chr. Monstein
"""

import datetime
import numpy as np
import time
import configparser
import time
import atexit

from utils import SunEphem, AntennaTracker, DISEQLogger
from utils import getlist
from utils import generate_times_for_range
from utils import generate_scanning_pattern
from utils import logger

#------------------------------------------------------------------------------
# Read configuration file
config = configparser.ConfigParser()
config.read('configsun.ini')
  
# Create Sun Ephem
sun_ephem = SunEphem(config.get('Location','Longitude'),
                    config.get('Location','Latitude'),
                    config.getfloat('Location','Elevation'),
                    config.getfloat('Location','Temperature'),
                    config.getfloat('Location','Pressure'))

#------------------------------------------------------------------------------
# Tracker parameter
maxrange  = config.getfloat('Tracker','MaxRange')
azimin    = config.getfloat('Tracker','AziMin')
azimax    = config.getfloat('Tracker','AziMax')
elemin    = config.getfloat('Tracker','EleMin')
elemax    = config.getfloat('Tracker','EleMax')
azipark   = config.getfloat('Tracker','AziPark')
elepark   = config.getfloat('Tracker','ElePark')
aziref    = config.getfloat('Tracker','AziRef')
eleref    = config.getfloat('Tracker','EleRef')
azidir    = config.getfloat('Tracker','AziDir')
eledir    = config.getfloat('Tracker','EleDir')
planecorr = config.getfloat('Tracker','PlaneCorr')
comport = config.get('Tracker','MyComport')

#------------------------------------------------------------------------------
# Calibration parameter
docali = config.getboolean('Calibration','DoCali')
azicold  = config.getfloat('Calibration','AziCold')
elecold  = config.getfloat('Calibration','EleCold')
azihot   = config.getfloat('Calibration','AziHot')
elehot   = config.getfloat('Calibration','EleHot')
try:
    timescali = config.get('Calibration','TimesCali')
except configparser.NoOptionError:
    tstartcali = config.get('Calibration','TStartCali')
    intervalcali = config.get('Calibration','IntervalCali')
    tstopcali = config.get('Calibration','TStopCali', fallback = "23:59")
    timescali = generate_times_for_range(tstartcali, tstopcali, intervalcali)
    
#------------------------------------------------------------------------------
# Scanning parameter
doscanning = config.getboolean('Scanning','DoScanning')
dazi = np.array(getlist(config.get('Scanning','Dazi')), dtype=np.float32)
dele = np.array(getlist(config.get('Scanning','Dele')), dtype=np.float32)
try:
    timesscanning = config.get('Scanning','TimeScanning')
except configparser.NoOptionError:
    tstartscanning = config.get('Scanning','TStartScanning')
    intervalscanning = config.get('Scanning','IntervalScanning')
    tstopscanning = config.get('Scanning','TStopScanning', fallback = "23:59")
    timesscanning = generate_times_for_range(tstartscanning, tstopscanning, intervalscanning)

scanning_pattern = generate_scanning_pattern(dazi, dele)

#------------------------------------------------------------------------------
# Create logger
dat = datetime.datetime.now().strftime("%Y-%m-%d") # "%Y-%m-%d %H:%M:%S"
log_filename ='DISEQ-' + dat + '-Sun.txt'
diseqlogger = DISEQLogger(log_filename, aziref, eleref)

#------------------------------------------------------------------------------
# Create antenna tracker
tracker = AntennaTracker(comport, 
                         azidir, aziref,
                         eledir, eleref,
                         planecorr)
atexit.register(tracker.close)  # Ensure port is closed on exit
        
# Initialize
isParked = False
did_last = None
while True: # main infinite loop
    # Get time in UTC
    dt = datetime.datetime.now(datetime.timezone.utc) 
    
    # Get azi, ele of the sun
    azisun, elesun = sun_ephem.get_sun_angles(dt)
    
    # Rotorangles that correspond to sun position
    rotorazi, rotorele = tracker.rotor_angles(azisun, elesun)
    
    # Now check what to do
    # During day: either calibration or scanning and if 
    # none of the two sun tracking
    # During night: parking
    
    # Check if parking
    if ((azisun > azimax) or (azisun < azimin) or 
        (elesun < elemin) or (elesun > elemax) or
        (np.abs(rotorele) > maxrange)    or
        (np.abs(rotorazi) > maxrange)):
        if (isParked == False):
            logger.info(f"Parking at Azimuth  = {azipark:7.2f}  Elevation   = {elepark:7.2f}")
            #tracker.send_position(azipark, elepark, maxrange)
            diseqlogger.log_execution(azipark, elepark,"Sun position possibly out of range")
            isParked = True
    else:
        isParked = False
        # Now check if calibration or scanning
        dt_time = dt.strftime("%H:%M")
        
        if dt_time in timescali and did_last != "calib": # Calibration
            logger.info(f"{dt.strftime("%Y-%m-%d %H:%M")}: performing calibration")
            diseqlogger.log_execution(azicold, elecold,'Cold sky reference position (Tref = Tcold)')
            #tracker.send_position(azicold, elecold)
            time.sleep(60)
            diseqlogger.log_execution(azicold, elecold,'Hot reference position e.g. on ground (To=Thot)')
            #tracker.send_position(azihot, elehot)
            time.sleep(60)
            did_last = "cali"
        
        elif dt_time in timesscanning and did_last != "scanning": # Scanning
            logger.info(f"Sun data: Azimuth   = {azisun:7.2f}  Elevation   = {elesun:6.2f}")
            logger.info(f"{dt.strftime("%Y-%m-%d %H:%M")}: performing scanning")
            for azi_offset, ele_offset in scanning_pattern:
                diseqlogger.log_execution(azisun + azi_offset, elesun + ele_offset, 'New scanning position')
                #tracker.send_position(azisun + azi_offset, elesun + ele_offset)
                time.sleep(3)
            did_last = "scanning"
        else: # Normal sun tracking
            #tracker.send_position(azisun, elesun)
            if did_last != "tracking":
                logger.info(f"Sun data: Azimuth   = {azisun:7.2f}  Elevation   = {elesun:6.2f}")
                logger.info(f"{dt.strftime("%Y-%m-%d %H:%M")}: performing tracking")
                diseqlogger.log_execution(azisun, elesun,'Pointing at the Sun')
            did_last = "tracking"
    time.sleep(10)