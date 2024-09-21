#%%[markdown]
"""
# Scope
This is a file for modelling and processing a solar tower
"""

#%%
import pathlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline, BSpline

from CombiCSP.misc import  heatmap2d
from CombiCSP import SolarTroughCalcs, SolarTowerCalcs, HOYS_DEFAULT,SolarSystemLocation
# import CombiCSP.SolarGeometry as sgh
from CombiCSP.storage import Tr

# from Demand_supply import *
#%% Load meteorological data from file and constants
hoy = HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
df_pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = df_pvgis.loc[:,'Gb(n)']

#%% Set Site location
sslCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt_hr=+2, alt=0)

#%% Tower related dimensions
# Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
# alt = 200*10e-3 #Height above sea level [m] # TODO this is probably 200*1e-3
# Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
# A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project

stc =  SolarTowerCalcs(alt_ = 200*10e-3 , Ht_km = 0.1, 
        Ar_m2 = 99.3 , A_helio_m2 = 225000,
        slobj=sslCrete)
oTow = stc.perform_calc(Ib)
#%%
plt.plot(hoy, oTow.data, label='1')
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()


#%%

oTow._df
# %%
a = oTow._df
# Conver Hour of Year to Hour of Day
a['hod']  = a.HOY % 24
a.head(30)
plt.plot(a.hod, a.Power_MW, '.')

# %%
# plot a bivariate kde for x hour of data, y power
import seaborn as sns
sns.kdeplot(x=a.hod, y=a.Power_MW, cmap='viridis', fill=True, cbar=True)

# %%
plt.plot(a.Ib_n, a.Power_MW, '.')
plt.xlabel('Irradiance (W/m2)') 
plt.ylabel('Power (MW)')
plt.title('Tower Power vs Irradiance')
# %%
# create a subset with Ib_n > 0
a_g0 = a[a.Ib_n > 0]
sns.kdeplot(x=a_g0.Ib_n, y=a_g0.Power_MW, cmap='viridis', fill=True, cbar=True)
plt.xlabel('Irradiance ($W/m^2$)')
plt.ylabel('Power (MW)')
plt.grid()
# %%
