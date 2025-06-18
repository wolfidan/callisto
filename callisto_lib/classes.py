import ephem
import math
import datetime
import os
import serial
import time
import logging

# Local imports 
import constants

class PythonLogger(logging.Logger):
    def __init__(self, processname, loglevel = "info"):
        # Initialize the parent Logger with a name
        super().__init__(name=processname, 
                         level=getattr(logging, loglevel.upper()))

        now = datetime.datetime.strftime(datetime.datetime.now(), '%Y%m%d')
        
        path_logs = os.path.join(constants.PATH_LOG, processname)
        if not os.path.exists(path_logs):
            os.makedirs(path_logs)
            
        logfile = os.path.join(path_logs,  f'{processname}_{now}.log')

        # Create a file handler
        handler = logging.FileHandler(logfile)
        handler.setLevel(getattr(logging, loglevel.upper()))

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)

        self.addHandler(handler)
        
class DISEQLogger(object):
    #------------------------------------------------------------------------------
    # Logging events to DISEQ file
    def __init__(self, filename, aziref, eleref):
        self.filename = filename
        if not os.path.exists(filename):
            with open(filename, "a") as fp:   # Save x/y-data in file
                fp.write('DiSEqC radio telescope tracker\n') # write header information
                fp.write('Author: Christian Monstein, HB9SCT\n') # write header information
                fp.write('Version: 2023-04-27\n') # write header information
                fp.write('\n')
                hms = datetime.datetime.now().strftime("%H:%M:%S") # "%Y-%m-%d %H:%M:%S"
                Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600
                fp.write('TIME [UT], Azimuth[deg], Elevation[deg], Message\n') # write header information
                st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(aziref) + '{:8.3f}, '.format(eleref) + 'Reference positions.\n'
                fp.write(st + '\n') 
    def log_execution(self, x,y, msg):
        with open(self.filename, "a") as fp:   # Save x/y-data in file
            hms = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3] # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600 + float(hms[9:])/3600000
            st = '{:8.5f},'.format(Tdec) + '{:8.3f},'.format(x) + '{:8.3f}, '.format(y) + msg
            fp.write(st+'\n') 
                   
class SunEphem(object):
    """
    Compute sun position and related angles using PyEphem.

    Parameters
    ----------
    longitude : float
        Longitude of the observation point in degrees (positive east).
    latitude : float
        Latitude of the observation point in degrees (positive north).
    elevation : float
        Elevation of the observation point in meters above sea level.
    pressure : float, optional
        Atmospheric pressure in millibars. If not provided, PyEphem uses its default.
    temperature : float, optional
        Ambient temperature in degrees Celsius. If not provided, PyEphem uses its default.

    Attributes
    ----------
    obs : ephem.Observer
        PyEphem observer instance configured with the given parameters.
    sun : ephem.Sun
        PyEphem Sun object for computing solar positions.
    """

    def __init__(self, longitude, latitude, elevation, 
                 pressure=None, temperature=None):
        self.obs = ephem.Observer()
        self.obs.lon = str(longitude)
        self.obs.lat = str(latitude)
        self.obs.elev = float(elevation)
        if temperature:
            self.obs.temp = float(temperature)
        if pressure:
            self.obs.pressure = float(pressure)
        self.sun = ephem.Sun()

    def get_sun_angles(self, dt):
        """
        Compute sun azimuth, elevation, hour angle, and declination for a given datetime.

        Parameters
        ----------
        dt : datetime.datetime
            Datetime for which to compute the sun angles.

        Returns
        -------
        azi : float
            Azimuth angle of the sun in degrees (0° = North, increasing clockwise).
        ele : float
            Elevation angle of the sun in degrees (0° = horizon, +90° = zenith).
        """
        self.obs.date = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}'.format(
            dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
        self.sun.compute(self.obs)
        azi = math.degrees(self.sun.az)
        ele = math.degrees(self.sun.alt)
    
        return azi, ele


class AntennaTracker(object):
    def __init__(self, logger, comport, 
                 azidir, aziref, eledir, eleref, 
                 planecorr = 0):
        
        self.logger = logger
        self.azidir = azidir
        self.eledir = eledir
        self.aziref = aziref
        self.eleref = eleref
        self.planecorr = planecorr
        
        # Open antenna tracker
        try:
            self.DiSEqC = serial.Serial(
                        port     = comport,
                        baudrate = 9600,
                        bytesize = serial.EIGHTBITS,
                        parity   = serial.PARITY_NONE,
                        timeout  = 2)
        except serial.SerialException as e:
            self.logger.error("Problem communication with tracker!")
            self.logger.error(e)
            raise
        
    def rotor_angles(self, azi, ele):
        rotorazi = self.azidir*(azi - self.aziref) # conversion 0° ... 360° -> +/- 180°
        rotorele = self.eledir*(ele - self.eleref) # conversion 0° ... 90° -> +/- 45°
        rotorele = (rotorele - self.planecorr*rotorazi)  # +/-1.25°
        return rotorazi, rotorele
    
    def send_position(self, azi, ele, maxrange = None):
        rotorazi, rotorele = self.rotor_angles(azi, ele)
        # Send a new position to antenna tracker
        if maxrange:
            cmd = 'max{:6.2f}\r'.format(maxrange) 
            self.DiSEqC.write(cmd.encode()) # set MaxRange
            time.sleep(1)
            
        cmd = 'azi{:8.3f}\r'.format(rotorazi) 
        self.DiSEqC.write(cmd.encode()) # set azimuth or hour angle
        time.sleep(1)
        
        cmd = 'ele{:8.3f}\r'.format(rotorele) 
        self.DiSEqC.write(cmd.encode()) # set elevation or declination
        time.sleep(1)
        
    def close(self):
        """
        Close the serial connection cleanly.
        """
        if hasattr(self, "DiSEqC") and self.DiSEqC.is_open:
            try:
                self.DiSEqC.close()
                self.logger.info("Serial connection closed.")
            except serial.SerialException as e:
                self.logger.warning(f"Failed to close serial port: {e}")
                