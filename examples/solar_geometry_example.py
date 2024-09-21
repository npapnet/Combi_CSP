#%% [markdown]
"""
# Solar Geometry Example
This is a script to investigate the solar geometry functions provided by the code

This was done to help understand better the current solar geometry function implementation 
and their results and to help debug them properly. 

**VALIDATED** Tests are currently required, because 
the currently available ones are the results from the initial implementaition  (sort of bootstrapping during 
the refactoring of the code from procedural to OOP). The following functions require validated tests:

- tsol (solar time)
- W (solar hour angle)
- z (zenith angle)
- ele (elevation angle)
- azim (azimuth angle)

TODO: complete test_SolarGeometry.py with the above functions and validate them.

"""

#%%
import pathlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline, BSpline

from CombiCSP import HOYS_DEFAULT,SolarSystemLocation

#%% Load data and constants
hoy = HOYS_DEFAULT

#%% Set Site location
sslCrete = SolarSystemLocation(lat=35, lon=24,  dt_gmt_hr=+2)

#%%
hoy2d =  np.arange(0,48, step=1)
solar_time = sslCrete.tsol(hoy=hoy2d)
plt.plot(hoy2d, solar_time-hoy2d)
plt.xlabel('Hour of the year')
plt.ylabel('(Solar Time - HOY) [hrs]')
#%%
solar_time