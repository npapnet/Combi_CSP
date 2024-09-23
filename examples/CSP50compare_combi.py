# %% [markdown]
"""---
title: CP50 example
author: N. Papadakis
date: 21/Sep/24
numbersections: true
---
"""
#%%[markdown]
"""
# Example of comparing the economics of a combined system ( solar tower and a solar trough system)

This one uses the solar_combi_system::CSPSystemCombination
to perform the calculation for a combination system

"""
#%%
import pathlib
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import numpy_financial as npf

from CombiCSP import HOYS_DEFAULT, SolarSystemLocation, SolarTowerCalcs, OutputContainer, Economic_environment, SolarTroughCalcs
from CombiCSP import CSPSystemCombination

import CombiCSP.misc as cspm
import CombiCSP.economics as cspe

from CombiCSP.storage import Tr

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

stc =  SolarTowerCalcs(alt_ = 200*10e-3 , Ht_km = 0.1, 
        Ar_m2 = 99.3 , A_helio_m2 = 225000,
        slobj=sslCrete)
# oTow = stc.perform_calc(Ib,transmittance=1)

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
#%% [markdown]
#  ## Define combined system
#%%

scmb = CSPSystemCombination(tow_lst=[stc], trough_lst=[strc])

#%% Perform analysis 
tow_args = {'transmittance':1, 'nG':0.97}
trough_args= {'Tr': 318, 'alignment':"NS"}
scmb.perform_calc(hoy=HOYS_DEFAULT,Ib=Ib, tow_args=tow_args, trough_args=trough_args)
#%%
#%%

#%% [markdown] ====================================================================================
# The following uses a class to perform the analysis. 
#%%

ee = Economic_environment(
    oil_price=oil_price, Eoil=Eoil,
    currency_units='USD'
)
#%%
area_list3 = []
cash_flow_list3 = []
tow_scenaria3 = []
for A_helio in np.arange(75000,125001,10000): # 100MW np.arange(150000,250001,10000):
    o_tmp = stc.mutate(A_helio=A_helio).perform_calc(Ib)
    tmp_res_Dic =ee.economics_for_Solar_tower(
            oTow= o_tmp,
            csp_area_costs= csp_area_costs,
            csp_energy_price=csp_energy_price,
            csp_discount_rate= csp_discount_rate,
            power_block_cost=power_block_cost,
        lifetime=range(30))
    area_list3.append(tmp_res_Dic['scenario_params']['A_helio_m2'] )
    cash_flow_list3.append(tmp_res_Dic['cash_flow'])
    tow_scenaria3.append(tmp_res_Dic['scenario_financial'])




ns_area_list = []  
ns_cash_flow_list = []
ns_trough_scenaria = []
for N_i in np.arange(800,1301,100): # 100MW np.arange(1000,2001,100):
    oTr = strc.mutate(N=N_i).perform_calcs_NS(Ib=Ib,hoy=HOYS_DEFAULT, Tr=Tr)
    tmp_res_Dic =ee.economics_for_SolarTrough(
            oTr= oTr,
            csp_area_costs= csp_area_costs,
            csp_energy_price=csp_energy_price,
            csp_discount_rate= csp_discount_rate,
            power_block_cost=power_block_cost,
        lifetime=range(30))
    ns_area_list.append(tmp_res_Dic['scenario_params']['area_m2'] )
    ns_cash_flow_list.append(tmp_res_Dic['cash_flow'])
    ns_trough_scenaria.append(tmp_res_Dic['scenario_financial'])

#%%

#%%

strc =  SolarTroughCalcs(foc_len=foc_len, N=1800, 
        L=L, Wr=Wr, Wc=Wc, Ws = Ws,
        slobj=sslCrete)
oTr = strc.perform_calcs_NS(Ib=Ib,hoy= HOYS_DEFAULT, Tr=Tr)
tmp_res_Dic =ee.economics_for_SolarTrough(
        oTr= oTr,
        csp_area_costs= csp_area_costs,
        csp_energy_price=csp_energy_price,
        csp_discount_rate= csp_discount_rate,
        power_block_cost=power_block_cost,
    lifetime=range(30))
#%%

# DPB = []
# for (x,y) in zip(area_list,cash_flow_list):
#     DPB.append(cspe.discounted_payback_period(csp_discount_rate, y).round(2))
#     #title(z)
#     #show()
# dpb_table = pd.DataFrame(DPB, index=area_list, columns = ['DPB (years)'])
# dpb_table

# T+NS optimum
A_helio_optNS = 125000
N_opt_NS = 800

stc_opt =  SolarTowerCalcs(alt_ = 200*10e-3 , Ht_km = 0.1, 
        Ar_m2 = 99.3 , A_helio_m2 = A_helio_optNS,
        slobj=sslCrete)
oTow = stc_opt.perform_calc(Ib,transmittance=1)
strc_opt =  SolarTroughCalcs(foc_len=foc_len, N=N_opt_NS, 
        L=L, Wr=Wr, Wc=Wc, Ws = Ws,
        slobj=sslCrete)
oTr = strc_opt.perform_calcs_NS(Ib=Ib,hoy= HOYS_DEFAULT, Tr=Tr)
#%%
# pd.merge(oTow.data_df, oTr.data_df, how="inner", on='HOY', )
# merge oTow.data_df and oTr.data_df on 'HOY' column, and use as suffixes tow and tr
combiNS = pd.merge(oTow.data_df, oTr.data_df, how="inner", on='HOY', suffixes=('_tow', '_tr'))
combiNS['Power_MW'] = combiNS['Power_MW_tow'] + combiNS['Power_MW_tr']
combiNS.head()
#%%

combiNS_power_2d_np = np.vstack(combiNS['Power_MW'].values).reshape((365,24))
#%%
tmp_res_Dic =ee.economics_for_Combination(
        oTr= oTr,
        oTow = oTow,
        csp_area_costs= csp_area_costs,
        csp_energy_price=csp_energy_price,
        csp_discount_rate= csp_discount_rate,
        power_block_cost=power_block_cost,
        # capital_csp= capital_csp,
        lifetime=range(30))
#%% [markdown]
""" 
## East West trough orienation 
"""
#%%
# # T+EW optimum
# A_helio_optEW = 125000
# N_opt_EW = 800

oTrEW = strc_opt.perform_calcs_EW(Ib=Ib,hoy= HOYS_DEFAULT, Tr=Tr)
combiEW = pd.merge(oTow.data_df, oTr.data_df, how="inner", on='HOY', suffixes=('_tow', '_tr'))
combiEW['Power_MW'] = combiEW['Power_MW_tow'] + combiEW['Power_MW_tr']
combiEW_xyz_np = np.vstack(combiEW['Power_MW'].values).reshape((365,24))

tmp_res_Dic =ee.economics_for_Combination(
    oTr= oTrEW,
    oTow = oTow,
    csp_area_costs= csp_area_costs,
    csp_energy_price=csp_energy_price,
    csp_discount_rate= csp_discount_rate,
    power_block_cost=power_block_cost,
    # capital_csp= capital_csp,
    lifetime=range(30)
    )
# %%
tmp_res_Dic.keys()
# %%
print('=========== Financial scenario outputs =============')
for k,v in tmp_res_Dic['scenario_financial'].items():
    if isinstance(v, str):
        print(f"{k} : {v}")
    else:
        try:
            iter(v)
        except TypeError:
            print(f"{k} : {v}")
#%%
print('=========== Tower scenario params =============')
for k,v in tmp_res_Dic['scenario_params']['tower'].items():
    if isinstance(v, str):
            print(f"{k} : {v}")
    else:
        try:
            iter(v)
        except TypeError:
            print(f"{k} : {v}")
    
print('=========== Trough scenario params =============')
for k,v in tmp_res_Dic['scenario_params']['trough'].items():
    if isinstance(v, str):
            print(f"{k} : {v}")
    else:
        try:
            iter(v)
        except TypeError:
            print(f"{k} : {v}")
# %%
