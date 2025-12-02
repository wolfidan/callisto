import ephem
import math
import datetime
import os
import serial
import time
import logging
import sys

# Local imports 
from . import constants

class PythonLogger(logging.Logger):
    def __init__(self, processname=None, loglevel = "info"):
        if processname:
            # Initialize the parent Logger with a name
            super().__init__(name=processname, 
                            level=getattr(logging, loglevel.upper()))

            now = datetime.datetime.strftime(datetime.datetime.now(), '%Y%m%d')
            
            dir_logs = os.path.join(constants.DIR_LOG, processname)
            if not os.path.exists(dir_logs):
                os.makedirs(dir_logs)
                
            logfile = os.path.join(dir_logs,  f'{processname}_{now}.log')

            # Create a file handler
            handler = logging.FileHandler(logfile)
            handler.setLevel(getattr(logging, loglevel.upper()))

            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)

            self.addHandler(handler)
        else:
            # stdout logger
             # Create a logger that logs to stdout
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)

            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter(
                "%(asctime)s [%(levelname)s] %(message)s"
            ))

            self.addHandler(handler)

        
class DISEQLogger(object):
    #------------------------------------------------------------------------------
    # Logging events to DISEQ file
    def __init__(self, filename, aziref, eleref):
        dirname = os.path.dirname(filename)
        if not dirname: # Filename without directory
            dirname = os.path.join(constants.DIR_LOG, "DISEQ")
            filename = os.path.join(dirname, filename)
        os.makedirs(dirname, exist_ok = True)
        self.filename = filename
        if not os.path.exists(filename):
            with open(filename, "a") as fp:   # Save x/y-data in file
                fp.write('# DiSEqC radio telescope tracker\n') # write header information
                fp.write('# Author: Christian Monstein, HB9SCT\n') # write header information
                fp.write('# Version: 2023-04-27\n') # write header information
                fp.write(f'# Reference positions: Azimuth={aziref:8.3f} deg, Elevation={eleref:8.3f} deg\n') # write header information
                fp.write("\n")
                # header
                fp.write('TIME [UT];Real azimuth [deg];Real elevation [deg];Rotor azimuth [deg];Rotor elevation [deg]; Message\n')

    def log_execution(self, xreal, yreal, xrot, yrot, msg):
        with open(self.filename, "a") as fp:   # Save x/y-data in file
            hms = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3] # "%Y-%m-%d %H:%M:%S"
            Tdec = float(hms[0:2]) + float(hms[3:5])/60 + float(hms[6:8])/3600 + float(hms[9:])/3600000
            st = f"{Tdec:8.5f};{xreal:8.3f};{yreal:8.3f};{xrot:8.3f};{yrot:8.3f};{msg}"
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
    """
    Control interface for an antenna rotor with optional DiSEqC command logging.

    This class provides methods to convert sky azimuth/elevation angles into
    rotor-relative angles, send positioning commands to a DiSEqC-compatible
    antenna rotor, and optionally log all issued commands for later inspection.

    Parameters
    ----------
    comport : str
        Serial port name used to communicate with the DiSEqC rotor
        (e.g., ``"COM3"`` or ``"/dev/ttyUSB0"``).
    azidir : float
        Direction multiplier for azimuth rotation (+1 or -1 depending on rotor geometry).
    aziref : float
        Reference azimuth angle used for computing rotor-relative angles.
    eledir : float
        Direction multiplier for elevation rotation (+1 or -1).
    eleref : float
        Reference elevation angle used for computing rotor-relative angles.
    planecorr : float, optional
        Plane-correction factor applied to elevation based on azimuth movement.
        Defaults to ``0``.
    python_logger : object, optional
        Logger instance with ``info()``, ``warning()``, and ``error()`` methods.
        If ``None``, a default stdout logger is created.
    diseqlog_file : str or None, optional
        If provided, the path to a file used to log all DiSEqC movement commands.
        If ``None``, logging is disabled.

    Raises
    ------
    serial.SerialException
        If the serial port cannot be opened.
    """

    def __init__(self, comport,
                 azidir, aziref, eledir, eleref,
                 planecorr=0,
                 python_logger=None,
                 diseqlog_file=None):

        if not python_logger:
            # stdout logger
            python_logger = PythonLogger()

        self.logger = python_logger
        self.azidir = azidir
        self.eledir = eledir
        self.aziref = aziref
        self.eleref = eleref
        self.planecorr = planecorr

        # Prepare DISEqC log file (optional)
        self.diseqc_logger = None
        if diseqlog_file is not None:
            self.diseqc_logger = DISEQLogger(diseqlog_file, aziref, eleref)

        # Open DiSEqC serial interface
        try:
            self.DiSEqC = serial.Serial(
                port=comport,
                baudrate=9600,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                timeout=2
            )
        except serial.SerialException as e:
            self.logger.error("Problem communicating with tracker!")
            self.logger.error(e)
            raise

    def rotor_angles(self, azi, ele):
        """
        Convert sky azimuth/elevation angles into rotor-relative coordinates.

        Parameters
        ----------
        azi : float
            Sky azimuth angle in degrees.
        ele : float
            Sky elevation angle in degrees.

        Returns
        -------
        rotorazi : float
            Rotor-relative azimuth angle.
        rotorele : float
            Rotor-relative elevation angle with optional plane correction applied.
        """
        rotorazi = self.azidir * (azi - self.aziref)
        rotorele = self.eledir * (ele - self.eleref)
        rotorele = rotorele - self.planecorr * rotorazi
        return rotorazi, rotorele

    def send_position(self, azi, ele, maxrange=None, transition_time=1):
        """
        Send a targeting command to the antenna rotor and optionally log it.

        Parameters
        ----------
        azi : float
            Target sky azimuth angle in degrees.
        ele : float
            Target sky elevation angle in degrees.
        maxrange : float or None, optional
            Optional maximum allowed movement range. If provided, a ``maxNNN``
            DiSEqC command is issued before positioning.
        transition_time : float, optional
            Delay (in seconds) inserted after each command to allow mechanical
            settling. Default is ``1``.

        Returns
        -------
        rotorazi : float
            Rotor-relative azimuth angle actually commanded.
        rotorele : float
            Rotor-relative elevation angle actually commanded.

        Notes
        -----
        If a DiSEqC log file is configured, this method will append a log entry
        containing both sky and rotor-relative angles.
        """
        rotorazi, rotorele = self.rotor_angles(azi, ele)

        # optional range limit
        if maxrange is not None:
            cmd = f"max{maxrange:6.2f}\r"
            self.DiSEqC.write(cmd.encode())
            time.sleep(transition_time)

        # Azimuth
        cmd = f"azi{rotorazi:8.3f}\r"
        self.DiSEqC.write(cmd.encode())
        time.sleep(transition_time)

        # Elevation
        cmd = f"ele{rotorele:8.3f}\r"
        self.DiSEqC.write(cmd.encode())
        time.sleep(transition_time)

        # LOG movement
        if self.diseqc_logger:
            self.diseqc_logger.log_execution(
                azi,
                ele,
                rotorazi,
                rotorele,
                msg=f"Commanded position az={rotorazi:.3f}, el={rotorele:.3f}"
            )

        return rotorazi, rotorele

    def close(self):
        """
        Close the serial connection cleanly.

        Closes the underlying DiSEqC serial interface if it is currently open.

        Raises
        ------
        serial.SerialException
            If the port cannot be closed.
        """
        if hasattr(self, "DiSEqC") and self.DiSEqC.is_open:
            try:
                self.DiSEqC.close()
                self.logger.info("Serial connection closed.")
            except serial.SerialException as e:
                self.logger.warning(f"Failed to close serial port: {e}")
