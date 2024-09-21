
#%% [markdown]
'''
@author: N. Papadakis

# Scope

This is an example of comparison between the solar tower and solar trough technologies using an oop approach. 
This is a rewrite of the code in CSP50compar.py by G. Arnaoutakis which used a procedural approach.
'''

#%%
import pathlib
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import numpy_financial as npf

from CombiCSP import HOYS_DEFAULT, SolarSystemLocation, SolarTowerCalcs, OutputContainer, Economic_environment, SolarTroughCalcs
import CombiCSP.misc as cspm
import CombiCSP.economics as cspe

from CombiCSP.storage import Tr

from CSP50_common_econ import *   # this is a development hack to avoid duplication
#import pcm


#%%[markdown]
'''
# Code
## Set Site location

The follwoing sets the location of the site in Crete
'''
#%% 
sslCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt_hr=+2, alt=0)
sslCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt_hr=+2, alt=0)

hoy = HOYS_DEFAULT

#%%[markdown]
'''
## Read meteo data

REad the Irradiance data from a file (the code can download from pvgis, but this is for offline use and development purposes).
'''

#%%
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis_data = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') #Atherinolakos
Ib = pvgis_data.loc[:,'Gb(n)']
#Ib = ineichen().dni
capital_csp = 5000000

#%%[markdown]
'''
## Define the solar tower 

The following defines the solar tower and calculates the power output for the given irradiance data.

- Ar: receiver area [m2] pp.44 in Pacheco
- alt: Height above sea level [m]
- Ht: Tower height [km]


'''
#%%
# Tower dimensions
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m]
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
#R = 5

#%% Tower related dimensions
# Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
# alt = 200*10e-3 #Height above sea level [m] # TODO this is probably 200*1e-3
# Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
# A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project

stc =  SolarTowerCalcs(alt_ = 200*10e-3 , Ht_km = 0.1, 
        Ar_m2 = 99.3 , A_helio_m2 = 225000,
        slobj=sslCrete)
oTow = stc.perform_calc(Ib,transmittance=1)

#%% [markdown]
"""
## Economical Analysis of Solar Tower

The following line creates an instance of the Economic_environment class which contains the oil price and the currency units
and has methods that can be used to calculate the economic parameters of the solar tower.

The following lines iterate over a range of helistat areas and calculate the cash flow and the scenaria for each area.

The data are stores in lists for further analysis:
- area_list3: the heliostat area
- cash_flow_list3: the cash flow
- tow_scenaria3: the scenarios parameters

TODO: I need to change from tuple to dictionary for the output of the scenarios or create another container class.

"""
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
    area_list3.append(tmp_res_Dic['A_helio'] )
    cash_flow_list3.append(tmp_res_Dic['cash_flow'])
    tow_scenaria3.append(tmp_res_Dic['scenaria'])



#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Wr = 0.07 # tube outer diameter [m]
Wc = 5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
N = 800

strc =  SolarTroughCalcs(foc_len=foc_len, N=N, 
        L=L, Wr=Wr, Wc=Wc, Ws = Ws,
        slobj=sslCrete)

ns_area_list = []  
ns_cash_flow_list = []
ns_trough_scenaria = []
for N_i in np.arange(800,1301,100): # 100MW np.arange(1000,2001,100):
    oTr = strc.mutate(N=N_i).perform_calcs_NS(Ib=Ib,hoy= hoy, Tr=Tr)
    tmp_res_Dic =ee.economics_for_SolarTrough(
            oTr= oTr,
            csp_area_costs= csp_area_costs,
            csp_energy_price=csp_energy_price,
            csp_discount_rate= csp_discount_rate,
            power_block_cost=power_block_cost,
        lifetime=range(30))
    ns_area_list.append(tmp_res_Dic['A_helio'] )
    ns_cash_flow_list.append(tmp_res_Dic['cash_flow'])
    ns_trough_scenaria.append(tmp_res_Dic['scenaria'])

#%%

#%%

strc =  SolarTroughCalcs(foc_len=foc_len, N=1800, 
        L=L, Wr=Wr, Wc=Wc, Ws = Ws,
        slobj=sslCrete)
oTr = strc.perform_calcs_NS(Ib=Ib,hoy= hoy, Tr=Tr)
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
oTr = strc_opt.perform_calcs_NS(Ib=Ib,hoy= hoy, Tr=Tr)

combiNS_np = (oTow.data + oTr.data).values
combiNS_xyz_np = np.vstack(combiNS_np).reshape((365,24))

tmp_res_Dic =ee.economics_for_Combination(
        oTr= oTr,
        oTow = oTow,
        csp_area_costs= csp_area_costs,
        csp_energy_price=csp_energy_price,
        csp_discount_rate= csp_discount_rate,
        power_block_cost=power_block_cost,
        # capital_csp= capital_csp,
        lifetime=range(30))

# # T+EW optimum
# A_helio_optEW = 125000
# N_opt_EW = 800

oTrEW = strc_opt.perform_calcs_EW(Ib=Ib,hoy= hoy, Tr=Tr)
combiEW_np = (oTow.data + oTrEW.data).values
combiEW_xyz_np = np.vstack(combiEW_np).reshape((365,24))

tmp_res_Dic =ee.economics_for_Combination(
        oTr= oTrEW,
        oTow = oTow,
        csp_area_costs= csp_area_costs,
        csp_energy_price=csp_energy_price,
        csp_discount_rate= csp_discount_rate,
        power_block_cost=power_block_cost,
        # capital_csp= capital_csp,
        lifetime=range(30))
# # %%

# %%
