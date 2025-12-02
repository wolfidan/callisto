
import numpy as np
import re
import datetime
import os
import xarray as xr
import pandas as pd

from .io_utils import read_diseq_log, read_fit
from .io_utils  import find_fit_file_for_timestamp, decimal_ut_to_datetime, extract_measurement

def read_fit_diseq(diseq_path):
    """
    This high-level function reads a DISEQ log, and from its timestamps, gets the corresponding
    fit measurements from the relevant fit files.

    Depending whether the measurements can be represented on a grid (i.e.  total number of records = 
    number of unique elevations x number of unique azimuths), the output will be

    - if data is not gridded: a dictionary with keys "elevation", "azimuth", "timestamp", "frequency" and "data"
    - if data is gridded: an xarray dataset with dimensions "elevation","azimuth", "frequency" and variables "timestamp" and "power"

    For a given DISEQ timestamp t[i], the code will average all fit measurements that fall between t[i] and t[i+1]. The reference time in
    the output will be t[i].

    Parameters
    ----------
    diseq_path : str
        DISEQ log filename, expect to be in the format <path>/*YYYYmmddTHHMM*

    Returns
    -------
    dict or xarray.Dataset
    """

    df = read_diseq_log(diseq_path)
    if df.empty:
        return df
    
    # Guess date from filename
    basename = os.path.basename(diseq_path)
    m = re.search(r"(20\d{6}T\d{4})", basename)
    if not m:
        raise ValueError("Cannot determine date from DISEQ filename. Embed YYYYmmddTHHMM in it.")

    date = datetime.datetime.strptime(m.group(1), "%Y%m%dT%H%M").date()
    # Build absolute datetimes
    df["datetime_utc"] = df["TIME [UT]"].apply(lambda ut: decimal_ut_to_datetime(date, ut))

    # Determine which FIT file each row belongs to
    df["fit_file"] = df["datetime_utc"].apply(find_fit_file_for_timestamp)

    # Load FIT files once
    df["measurement"] = None

    for f in df["fit_file"].unique():
        if f is None:
            continue
        fit_data = read_fit(f)

        subset = df["fit_file"] == f
        ts0 = df.loc[subset, "datetime_utc"]
        ts1 = ts0.shift(-1)     # ts1[i] = ts0[i+1], last one becomes NaT
        values = [
            extract_measurement(fit_data, t0, t1)
            for t0, t1 in zip(ts0, ts1)
        ]
        df.loc[subset, "measurement"] = pd.Series(values, dtype="object").iloc[:].values

    nazi = len(np.unique(df["Real azimuth [deg]"]))
    nele = len(np.unique(df["Real elevation [deg]"]))
    frequencies = fit_data["freqs"]

    is_gridded = False
    if (nazi * nele) == len(df):
        is_gridded = True

    # Extract values
    az = df["Real azimuth [deg]"].to_numpy()
    az_unique = np.unique(az)
    el = df["Real elevation [deg]"].to_numpy()
    el_unique = np.unique(el)
    datetimes = df["datetime_utc"].to_numpy()

    # Measurement column: convert each list to np.array
    measurements = np.stack(df["measurement"].to_numpy())  # shape (Nrecords, Nfreq)
    if is_gridded:
        # Build xarray Dataset
        ds = xr.Dataset(
            {
                "power": (("azimuth", "elevation", "frequency"), 
                        measurements.reshape(nazi, nele, -1)),  # reshape needed
                "datetime": (("azimuth", "elevation"), 
                            datetimes.reshape(nazi, nele)),
            },
            coords={
                "azimuth": az_unique,
                "elevation": el_unique,
                "frequency": frequencies,
            }
        )
    else:
        ds = {"elevation": el, "azimuth": az, "datetime": datetimes, 
        "power": measurements, "frequency": frequencies}
    
    return ds
