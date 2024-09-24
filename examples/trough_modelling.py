#%%[markdown]
"""
---
title: Trough Modelling
author: N. Papadakis

---
# Trough modelling Script

This is an example for showcasing the functionality of the CombiCSP library with respect to the solar trough modelling.
"""
#%%
import pathlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy.interpolate import make_interp_spline, BSpline

from CombiCSP.misc import  heatmap2d
from CombiCSP import SolarTroughCalcs, HOYS_DEFAULT,SolarSystemLocation
# import CombiCSP.SolarGeometry as sgh
from CombiCSP.storage import Tr

# from Demand_supply import *
#%% Load data and constants
hoy = HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
df_pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = df_pvgis.loc[:,'Gb(n)']

#%% Set Site location
sslCrete = SolarSystemLocation(lat=35, lon=24, dt_gmt_hr=+2, alt=0)


#%%
# Trough dimensions
# foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
# Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
# L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
# N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
# Wr_ini=0.07 # tube outer diameter [m]
# Wc_ini=5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM


sotr = SolarTroughCalcs(
        foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
        ,N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
        ,L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  
        ,Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
        ,Wr = 0.07 # tube outer diameter [m]
        ,Wc = 5.76
        , slobj=sslCrete
        )
oew = sotr.perform_calc(Ib=Ib, Tr=Tr, alignment='EW')
ons = sotr.perform_calc(Ib=Ib, Tr=Tr, alignment='NS')
#%%
plt.plot(hoy, ons.data_df['Power_MW'], '.', alpha=0.3, label='NS')#,xlim(100,600)
plt.plot(hoy, oew.data_df['Power_MW'], '.', alpha=0.3,   label='EW' )#,xlim(100,600)
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.legend()
plt.title('Trough')
#xlim(0,87.60), ylim(0,80)
#np.savetxt('datah.txt',datah.T ,delimiter=',') #save transposed data
#np.savetxt('dataxyz.txt',datah.T ,delimiter=',') #save tro_xyz data
plt.show()

#%%

plt.title('Trough N-S')
heatmap2d(ons.data4surf().T)
plt.title('Trough E-W')
heatmap2d(oew.data4surf().T)


# # %%


# CombiCSP.misc.heatmap_sns(tro_xyz.T, title='Trough E-W')
# # %%
# %% [markdown]
"""
## Solar Trough EW Alignment

The following code block demonstrates the calculation of the solar trough power output for the East-West alignment.
"""
#%%
oew = sotr.perform_calc(Ib=Ib, Tr=Tr, alignment='EW')
#%
# %% [markdown]
"""
### Power output vs Hour of Day

The following code block demonstrates the power output of the solar trough for the East-West alignment vs the hour of the day.
"""
# %%
_a = oew.data_df
# Conver Hour of Year to Hour of Day
_a['hod']  = _a.HOY % 24
_a.head(30)
plt.plot(_a.hod, _a.Power_MW, '.')


# %% [markdown]
"""
### Solar Trough Power distribution vs Hour of Day

The following code block demonstrates the power output of a  solar Trough vs the hour of the day.
"""
# %%
# plot a bivariate kde for x hour of data, y power
sns.kdeplot(x=_a.hod, y=_a.Power_MW, cmap='viridis', fill=True, cbar=True)

# %% [markdown]
"""
### Solar Trough  Power distribution vs Irradiance

The following code block demonstrates the power output of a  solar tower vs the Irradiance: 

- The first plot is a scatter plot of the power output vs the irradiance. 
- The second plot is a bivariate kde plot of the power output vs the irradiance (for irradiance values greater than zero).
"""
# %%
plt.plot(_a.Ib_n, _a.Power_MW, '.')
plt.xlabel('Irradiance (W/m2)') 
plt.ylabel('Power (MW)')
plt.title('Solar Trough Power vs Irradiance')
# %%
# create a subset with Ib_n > 0
_a_g0 = _a[_a.Ib_n > 0]
sns.kdeplot(x=_a_g0.Ib_n, y=_a_g0.Power_MW, cmap='viridis', fill=True, cbar=True)
plt.xlabel('Irradiance ($W/m^2$)')
plt.ylabel('Power (MW)')
plt.grid()
# %%
#%% [markdown]
"""
# Financial assessment of Solar Trough investmet
The code is able to perform and report a quick financial assessment of the solar trough investment.

The financial assessment is based on the following inputs

- Capital investment cost of the solar trough (land and power block)
- sale price of electricity
- discount rate
- lifetime of the investment

The following are included indirectly  in the financial assessment
- Operating and maintenance costs (as 4% per of the capital investment per annum) 
"""
# %%
from CombiCSP import Economic_environment

oil_price = 60 # np.arange(12, 112, 10)# [$/barrel] https://www.statista.com/statistics/262860/uk-brent-crude-oil-price-changes-since-1976/
Eaux = 20 # [MWh]
csp_area_costs = 25+ 150+  60 # site dev, coll cost, htf cost
power_block_cost = 9.1e5 # [$/MW]
csp_energy_price = 248 # [$/MWh]
discount_rate = 0.09
Eoil = Eaux*0.5883 # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000


eenv = Economic_environment(
    oil_price=oil_price, 
    Eoil=Eoil,
    currency_units='USD'
)
#%%
stc_fin = sotr.financial_assessment(oTr=oew, 
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
stc_fin['cash_flow_df']
# %%
stc_fin['scenario_financial']
# %%
