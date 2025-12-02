import datetime
import glob
import os
from io import StringIO
import re

import xarray as xr
import numpy as np
import pandas as pd
import requests
from astropy.io import fits

# Local imports
from . import constants
from .time_utils import decimal_ut_to_datetime

def read_fit(path: str, cache: dict = {}) -> dict:
    """
    Read a FIT file and extract data, axes, metadata, and datetimes.

    Parameters
    ----------
    path : str
        Path to the FIT file.

    Returns
    -------
    dict
        Dictionary containing:
        - **data** : np.ndarray  
          Main data array.
        - **freqs** : np.ndarray  
          Frequency axis.
        - **timeax** : np.ndarray  
          Time axis (seconds since file start).
        - **dT** : float  
          Time resolution.
        - **date** : str  
          Observation date from header.
        - **T0** : str  
          Start time of observation.
        - **instrument** : str  
          Instrument name.
        - **content** : str  
          Content description.
        - **datetimes** : np.ndarray of datetime.datetime  
          Absolute datetimes with UTC timezone.
    """
    if path in cache:
        return cache[path]

    with fits.open(path) as hdu:
        d = {
            "data": hdu[0].data,
            "freqs": hdu[1].data["Frequency"][0],
            "timeax": hdu[1].data["Time"][0],
            "dT": hdu[0].header["CDELT1"],
            "date": hdu[0].header["DATE-OBS"],
            "T0": hdu[0].header["TIME-OBS"],
            "instrument": hdu[0].header["INSTRUME"],
            "content": hdu[0].header["CONTENT"],
        }

    # Build absolute datetime axis
    base = datetime.datetime.strptime(
              f"{d['date']}_{d['T0']}", "%Y/%m/%d_%H:%M:%S.%f"
    ).replace(tzinfo=datetime.UTC)

    d["datetimes"] = np.array(
        [base + datetime.timedelta(seconds=dt) for dt in d["timeax"]]
    )

    cache[path] = d
    return d



def read_diseq_log(filename: str, dir_log_motor: str | None = None) -> pd.DataFrame:
    """
    Read a motor DISEQ log

    Parameters
    ----------
    filename : str
        Filename (without directory)
    dir_log_motor : str or None, optional
        Base path to DISEQ logs. If None, uses `constants.DIR_LOG/DISEQ`.

    Returns
    -------
    pandas.DataFrame
        Parsed motor log.
    """
    if dir_log_motor is None:
        dir_log_motor = os.path.join(constants.DIR_LOG, "DISEQ")

    fullpath = os.path.join(dir_log_motor, filename)

    df = pd.read_csv(fullpath, sep=";", header=4)
    df.rename(columns=lambda x: x.strip(), inplace=True)
    df["Message"] = df["Message"].apply(lambda x: x.strip())
    return df


def get_meteoswiss_data(station: str) -> pd.DataFrame:
    """
    Download the OGDS-SMN *t_now* CSV file for the given station.

    Parameters
    ----------
    station : str
        Station name (case-insensitive).

    Returns
    -------
    pandas.DataFrame
        Station's *t_now* dataset.
    """
    station_l = station.lower()
    url = (
        f"https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn/"
        f"{station_l}/ogd-smn_{station_l}_t_now.csv"
    )

    resp = requests.get(url)
    resp.raise_for_status()

    return pd.read_csv(StringIO(resp.text), sep=";")


def get_goes_data() -> pd.DataFrame:
    """
    Fetch 7-day primary GOES X-ray flux data from NOAA.

    Returns
    -------
    pandas.DataFrame
        GOES X-ray dataset.
    """
    return pd.read_json(
        "https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json"
    )


def get_fit_files(day: str, path_fit_folder: str | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Get FIT files and their timestamps for a given day.

    Parameters
    ----------
    day : str
        Date string in `YYYYMMDD` format.
    path_fit_folder : str or None, optional
        Folder containing FIT files. If None, uses
        `constants.PATH_DATA/raw/FITfiles`.

    Returns
    -------
    tuple of (np.ndarray, np.ndarray)
        - **files** : ndarray of str  
          List of FIT filenames.
        - **times** : ndarray of datetime.datetime  
          Extracted timestamps (UTC).
    """
    if path_fit_folder is None:
        path_fit_folder = os.path.join(constants.PATH_DATA, "raw", "FITfiles")

    filename_pattern = "meteoswiss_{:s}_01.fit"
    search_path = os.path.join(path_fit_folder, filename_pattern.format(day + "_*"))

    files = np.array(sorted(glob.glob(search_path)))

    times = np.array(
        [
            datetime.datetime.strptime(
                os.path.basename(fn),
                filename_pattern.format("%Y%m%d_%H%M%S"),
            ).replace(tzinfo=datetime.UTC)
            for fn in files
        ]
    )

    return files, times


def getlist(option: str, sep: str = ",", chars: str | None = None) -> list[str]:
    """
    Split a ConfigParser-style option string into a list.

    Parameters
    ----------
    option : str
        Input string to split.
    sep : str, optional
        Separator, by default ",".
    chars : str or None, optional
        Characters to strip from each item.

    Returns
    -------
    list of str
        Cleaned string items.
    """
    return [chunk.strip(chars) for chunk in option.split(sep)]


def find_fit_file_for_timestamp(ts: datetime.datetime):
    """
    FIT files end with the END of the recording period:
        meteoswiss_YYYYMMDD_HHMMSS_01.fit

    We choose the first file whose timestamp >= ts.
    """
    flist = glob.glob(os.path.join(constants.DIR_DATA, "raw", "FITfiles", "meteoswiss_*.fit"))
    candidates = []

    for f in flist:
        base = os.path.basename(f)
        m = re.match(r"meteoswiss_(\d{8})_(\d{6})_", base)
        if not m:
            continue

        ds, ts6 = m.groups()
        start_time = datetime.datetime.strptime(ds + ts6, "%Y%m%d%H%M%S").replace(tzinfo=datetime.UTC)

        if start_time < ts:
            candidates.append((start_time, f))

    if not candidates:
        return None

    # Pick the one with the largest start_time
    candidates.sort(key=lambda x: x[0])
    return candidates[-1][1]


def extract_measurement(fit_data, ts0, ts1=None):
    """
    Extract measurement(s) from fit_data.
    
    Parameters
    ----------
    fit_data : dict or xarray-like
        Contains 'datetimes' (array of datetime objects)
        and 'data' of shape (freq, time).
    ts0 : pandas.Timestamp
        Start timestamp (UTC).
    ts1 : pandas.Timestamp, optional
        End timestamp (UTC). If None, returns measurement at
        the closest time to ts0. If provided, returns the 
        average of all measurements between ts0 and ts1.
        
    Returns
    -------
    ndarray
        If ts1 is None → measurement at closest time.
        If ts1 is provided → averaged measurement across range.
    """
    t = fit_data["datetimes"]
    t = np.asarray(t)

    # Convert timestamps to Python datetime for comparison
    ts0 = ts0.to_pydatetime()

    if not ts1 or pd.isna(ts1):
        # --- Default behavior: closest measurement to ts0 ---
        idx = np.argmin(np.abs(t - ts0))
        return fit_data["data"][:, idx]

    else:
        ts1 = ts1.to_pydatetime()
        # --- Average all measurements between ts0 and ts1 ---
        # Ensure correct ordering even if user gives ts1 < ts0
        t_start = min(ts0, ts1)
        t_end = max(ts0, ts1)

        # Boolean mask of times inside window
        mask = (t >= t_start) & (t <= t_end)

        # No data in interval → return NaNs
        if not np.any(mask):
            return np.full(fit_data["data"].shape[0], np.nan)

        # Average over time axis
        return fit_data["data"][:, mask].mean(axis=1)

