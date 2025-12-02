# -*- coding: utf-8 -*-
"""
This script does a raster scan of the sun, similarly to what the Leonardo DX50 weather radar is doing
Two operating modes are supported:
- coarse: scans the whole sky to find for the possible aziimuth and elevation of the sun
- fine: scans around the estimated aziimuth and elevation of the sun, to fine-tune the estimation of the coarse scan
"""

#%%
import numpy as np
from callisto_lib.readers import read_fit_diseq

ds = read_fit_diseq("DISEQ-20251202T1313-Sunraster-coarse.txt")

#%%
import matplotlib.pyplot as plt
target_freq = 10637
idx_frequency = np.where(ds["frequency"].astype(int) == target_freq)[0][0]
plt.pcolormesh(ds["azimuth"], ds["elevation"], ds["power"][:,:,idx_frequency]); plt.colorbar()