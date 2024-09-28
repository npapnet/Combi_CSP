#%%[markdown]

"""
---
title: Optimise LCOE for a CSP plant in Crete
author: N.Papadakis
---
# Scope

optimise the size of the plant based on different combinations of nominal power solar towers and solar trucks.

# Problem formulation

The problem is the following, assume nominal power of 50 MW and a solar radiation data for Crete, Ierapetra. 

The goal is to find the optimal combination of solar tower and solar troughs that will minimise the LCOE.

1. Find the number of units solar trough object for 50 MW (N_50)
2. Create a linspace from 0 to N_50 
3. For each value $N_i$ calculate the LCOE. The LCOE calculation has the following steps
    a. Create a solar trough object with $N_i$ untis
    b. Calcualate the nominal power outut
    c. Using the nominal power output calculate the LCOE, calculate the optimal area of the solar tower
    d. perform caluclations for the objects
    e. Calculate the total flows and cost 
    f. Calculate the LCOE

"""

#%%
import pathlib
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import numpy_financial as npf

from CombiCSP import HOYS_DEFAULT, SolarSystemLocation, SolarTowerCalcs, OutputContainer, EconomicEnvironment, SolarTroughCalcs
from CombiCSP import CSPSystemCombination

import CombiCSP.misc as cspm
import CombiCSP.financial.economics as cspe
from CombiCSP.financial.lcoe import LCOECalculator

from CSP50_common_econ import *   # this is a development hack to avoid duplication

#%% [markdown]
"""
## Set the site
Set the Site location 
at Crete::Ierapetra(lat=35, lon=24
"""
#%% 
sslCrete = SolarSystemLocation(lat=35, lon=24,  dt_gmt_hr=+2, alt=0)

#%% [markdown]
# ## read the radiation datadata from local file
#%%
# 
# "tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis_data = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') #Atherinolakos
Ib = pvgis_data.loc[:,'Gb(n)']
#Ib = ineichen().dni
# capital_csp = 5000000


#%%[markdown]
# ## Define trough object
# Dimensions that are not changing:
# - foc_len : [m] focal length CSPP T.1 in Mosleh19
# - Wr : tube outer diameter [m]
# - Wc: collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
# - Ws : [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
# - L : [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
#
# The optimisatoin variable is:
# - N = 800 # Number of units this is the 
#
#  Define solar trough object **strc**

#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Wr = 0.07 # tube outer diameter [m]
Wc = 5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
N = 800 # Number of units this is the optimisation variable

strc =  SolarTroughCalcs(foc_len=foc_len, N=N, 
        L=L, Wr=Wr, Wc=Wc, Ws = Ws,
        slobj=sslCrete)


#%%[markdown]
"""
## Define tower object
Dimensions that are not changing:

- receiver area [$m^2$]
- altitude [km]
- Tower Height km]
"""
# Define solar tower object oTow
#%%
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m] #FIXME this si probably [km] instead of [m]
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
# A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project
# Define object Tower related dimensions **stc**

stoc =  SolarTowerCalcs(alt_ = 200*10e-3 , Ht_km = 0.1, 
        Ar_m2 = 99.3 , A_helio_m2 = 225000,
        slobj=sslCrete)
# oTow = stc.perform_calc(Ib,transmittance=1)

#%%[markdown]
"""
create a function that will calculate the LCOE for a given number of units of the solar trough object
"""
#%%
def process_one_scenario(no_trough_units_i = 1000, total_power = 50, verbose:bool=False):
    
    T_amb_C = 15
    nG= 0.97

    Eaux = 20 # [MWh]
    oil_price = 60 # np.arange(12, 112, 10)# [$/barrel] https://www.statista.com/statistics/262860/uk-brent-crude-oil-price-changes-since-1976/
    Eoil = Eaux*0.5883 # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000
    Eaux = 20 # [MWh]
    csp_area_costs = 25+ 150+  60 # site dev, coll cost, htf cost
    power_block_cost = 9.1e5 # [$/MW]
    csp_energy_price = 248 # [$/MWh]
    discount_rate = 0.09

    
    ee = EconomicEnvironment(
        oil_price=oil_price, Eoil=Eoil,
        currency_units='USD'
    )

    # find the nominal power of the solar trough with N number of units
    strc_i = strc.mutate(N=no_trough_units_i)
    nominal_trough_power = strc_i.calculate_nominal_power_MW(Tr=318, T_amb=T_amb_C, nG=0.97)

    # find the optimal area of the solar tower
    target_tower_power_MW = total_power - nominal_trough_power
    heliostat_area_i = stoc.find_area_for_nominal_power(target_power_MW=target_tower_power_MW, T_r_C=565)
    stoc_i = stoc.mutate(A_helio=heliostat_area_i)

    # check total power
    toc_i_power =  stoc_i.calculate_nominal_power_MW(T_r_C=565, T_amb_C=T_amb_C, nG=0.97)
    tro_i_power = strc_i.calculate_nominal_power_MW(Tr=318, T_amb=T_amb_C, nG=0.97) 
    total_power_i = tro_i_power + toc_i_power
    assert np.isclose(total_power_i, total_power, atol=1e-3, rtol=1e-4), f"Total power {total_power_i} MW, while expecting {total_power} MW"	
    print(f"Total power {total_power_i} MW")

    # strc_i.perform_calc(Ib,transmittance=1)
    oTow_i = stoc_i.perform_calc(Ib,transmittance=1)
    oTro_i = strc_i.perform_calc(Ib=Ib, Tr=318, alignment="NS")

    # calculate total power 
    total_power =  oTow_i.data_df['Power_MW'] + oTro_i.data_df['Power_MW']

    if verbose:
        print(f" Tower:  Nominal power  ={stoc_i.calculate_nominal_power_MW(T_r_C=318, T_amb_C=T_amb_C, nG=0.97):.6g} MW")
        print(f" Tower:  Max     power  ={oTow_i.PowerMax_MW:.6g} MW")

        print(f" Trough:  Nominal power  ={strc_i.calculate_nominal_power_MW(Tr=318, T_amb=15,nG=0.97):.6g} MW")
        print(f" Trough:  Max     power  ={oTro_i.PowerMax_MW:.6g} MW")


    tow_capical_cost = stoc_i.A_helio_m2* csp_area_costs + stoc_i.calculate_nominal_power_MW()*power_block_cost
    tow_opex = 0.04*tow_capical_cost

    tro_capical_cost = strc_i.area * csp_area_costs + strc_i.calculate_nominal_power_MW()*power_block_cost
    tro_opex = 0.04*tro_capical_cost

    if verbose:
        print(f" Tower  :  Capital cost = {tow_capical_cost:.6g} USD")
        print(f"        :  OPEX         = {tow_opex:.6g} USD")
        print(f" Trough :  Capital cost = {tro_capical_cost:.6g} USD")
        print(f"        :  OPEX         = {tro_opex:.6g} USD")

    # Create df with costs and energy
    year_id = np.arange(0,31)
    costs_arr = np.ones(shape=(31,))* (tow_opex+tro_opex)
    costs_arr[0] = (tow_capical_cost + tro_capical_cost)
    energy = np.ones(shape=(31,))* ( oTow_i.Energy_MWh + oTro_i.Energy_MWh)

    df = pd.DataFrame(data={'year':year_id, 'costs_per_annum':costs_arr, 'energy_MWh':energy})
    df.head()

    lcoe_calc = LCOECalculator()
    lcoe_val = lcoe_calc.lcoe_Eu_MWh(df, discount_rate)

    results_dic = {
        "lcoe USD/MWh": lcoe_val,
        "tower nominal power MW": stoc_i.calculate_nominal_power_MW(T_r_C=565, T_amb_C=T_amb_C, nG=0.97),
        "tower max power MW": oTow_i.PowerMax_MW,
        "tower energy MWh": oTow_i.Energy_MWh,
        "trough nominal power MW": strc_i.calculate_nominal_power_MW(Tr=318, T_amb=15,nG=0.97),
        "trough max power MW": oTro_i.PowerMax_MW,
        "trough energy MWh": oTro_i.Energy_MWh,
    }
    return results_dic

#%%
total_power = 500
no_units_50 = strc.find_no_units_for_nominal_power_MW(target_MW=total_power ,Tr=318)
print(f"Number of units for 50 MW is {no_units_50}")
res_list = []
for i in range(1, int(no_units_50), 100):
    print(f"Processing scenario {i}")
    results = process_one_scenario(i, total_power=total_power )
    res_list.append(results)
    print(results)

df_results = pd.DataFrame(res_list)
df_results.head()
df_results.describe()
#%%
import matplotlib.pyplot as plt 
plt.plot( df_results['tower nominal power MW'], df_results['lcoe USD/MWh'], '.',label='LCOE')
plt.plot( df_results['tower nominal power MW'].iloc[[0,-1]] , df_results['lcoe USD/MWh'].iloc[[0,-1]], 'r--', alpha=0.5, label='straight line' )
plt.xlabel('Tower nominal power MW')
plt.ylabel('LCOE USD/MWh')

#%%
