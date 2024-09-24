#%%[markdown]
"""
# Scope
This is a file for modelling and processing a solar tower
"""
#%%
%load_ext autoreload
%autoreload 2

import pathlib
import numpy as np
import pandas as pd


from scipy.interpolate import make_interp_spline, BSpline
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

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
sslCrete = SolarSystemLocation(lat=35, lon=24,dt_gmt_hr=+2, alt=0)

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
plt.plot(hoy, oTow.data_df, label='1')
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()


#%%

oTow.data_df
# %% [markdown]
"""
### Solar Tower Power output vs Hour of Day

The following code block demonstrates the power output of the solar trough for the East-West alignment vs the hour of the day.
"""
# %%
a = oTow.data_df
# Conver Hour of Year to Hour of Day
a['hod']  = a.HOY % 24
a.head(30)
plt.plot(a.hod, a.Power_MW, '.')

# %% [markdown]
"""
### Solar Tower Power distribution vs Hour of Day

The following code block demonstrates the power output of a  solar tower vs the hour of the day.
"""
# %%
# plot a bivariate kde for x hour of data, y power
sns.kdeplot(x=a.hod, y=a.Power_MW, cmap='viridis', fill=True, cbar=True)

# %% [markdown]
"""
### Solar Tower Power distribution vs Hour of Day

The following code block demonstrates the power output of a  solar tower vs the hour of the day.
"""
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
oTow.data_df
#%% [markdown]
"""
# Financial assessment of Solar investmet
The code is able to perform and report a quick financial assessment of the solar tower investment.

The financial assessment is based on the following inputs

- Capital investment cost of the solar tower (land and power block)
- sale price of electricity
- discount rate
- lifetime of the investment

The following are not included in the financial assessment
- Operating and maintenance costs
"""
# %%
from CombiCSP import EconomicEnvironment

oil_price = 60 # np.arange(12, 112, 10)# [$/barrel] https://www.statista.com/statistics/262860/uk-brent-crude-oil-price-changes-since-1976/
Eaux = 20 # [MWh]
csp_area_costs = 25+ 150+  60 # site dev, coll cost, htf cost
power_block_cost = 9.1e5 # [$/MW]
csp_energy_price = 248 # [$/MWh]
discount_rate = 0.09
Eoil = Eaux*0.5883 # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000


eenv = EconomicEnvironment(
    oil_price=oil_price, 
    Eoil=Eoil,
    currency_units='USD'
)
#%%
df = stc_fin = stc.financial_assessment(oTow= oTow, 
        ee=eenv,
        csp_area_costs= csp_area_costs,
        power_block_cost= power_block_cost,
        csp_energy_price= csp_energy_price,
        csp_discount_rate= discount_rate,
        lifetime=30)
# %%
stc_fin
# %%
stc_fin.keys()
# %%
stc_fin['cash_flow']
# %%
stc_fin['scenario_financial']
# %%
