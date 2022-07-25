# new OOP example of comparison
#%%
import pathlib
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import numpy_financial as npf

from CombiCSP import HOYS_DEFAULT, SolarSystemLocation, SolarTowerCalcs, OutputContainer, Economic_environment, SolarTroughCalcs
import CombiCSP.SolarGeometry as sgh
import CombiCSP.misc as cspm
import CombiCSP.economics as cspe

from CombiCSP.storage import Tr

from CSP50_common_econ import *   # this is a development hack to avoid duplication
#import pcm
#%%

#%% Set Site location
sslCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)

hoy = HOYS_DEFAULT
#%% Imported from CSPEcon ==========================================================================================
# from CSPecon import csp_area_costs, power_block_cost, csp_energy_price \
#     ,csp_discount_rate, oil_price, Eoil

# from CSPCret import Ptrough, Etrough



# End of Import from CSPEcon ==========================================================================================
#%%
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis_data = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') #Atherinolakos
Ib = pvgis_data.loc[:,'Gb(n)']
#Ib = ineichen().dni
capital_csp = 5000000


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

stc =  SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, 
        Ar = 99.3 , A_helio = 225000,
        slobj=sslCrete)
oTow = stc.perform_calc(Ib,transmittance=1)

#%% [markdown]
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
            capital_csp=capital_csp,
        lifetime=range(30))
    area_list3.append(tmp_res_Dic['A_helio'] )
    cash_flow_list3.append(tmp_res_Dic['cash_flow_tow'])
    tow_scenaria3.append(tmp_res_Dic['tow_scenaria'])



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

# # T+NS optimum
# A_helio_optNS = 125000
# N_opt_NS = 800
# tower_opt = solarII(Ib,1,IAM_tow(hoy),A_helio_optNS,Ar)
# trough_opt = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_NS)
# combiNS = tower_opt + trough_opt
# combiNS_xyz = np.vstack(combiNS).reshape((365,24)) # reshape 8760,1 to 365,24
# plt.title('Tower + Trough N-S')
# cspm.heatmap2d(combiNS_xyz.T)

# area_combiNS = Ac(Wc, L, N_opt_NS)
# PcombiNS = np.amax(combiNS) # used in CSPecon .round(2)
# EcombiNS = integrate.trapz(combiNS).round(2) # used in CSPecon

# capital_combiNS = (A_helio_optNS+area_combiNS)*csp_area_costs + PcombiNS*power_block_cost
# revenue_combiNS = cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS)
# cash_flow_combiNS = [-capital_combiNS] + [revenue_combiNS for i in range(30)]
# dpb_combiNS = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiNS)
# npv_combiNS = npf.npv(csp_discount_rate, [-capital_combiNS] \
#     + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])
# irr_combiNS = npf.irr([-capital_csp] \
#     + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])



# # T+EW optimum
# A_helio_optEW = 125000
# N_opt_EW = 800
# tower = solarII(Ib,1,IAM_tow(hoy),A_helio_optEW,Ar)
# troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_EW)
# combiEW = tower + troughew
# combiEW_xyz = np.vstack(combiEW).reshape((365,24)) # reshape 8760,1 to 365,24
# plt.title('Tower + Trough E-W')
# cspm.heatmap2d(combiEW_xyz.T)

# area_combiEW = Ac(Wc, L, N_opt_EW)
# PcombiEW = np.amax(combiEW) # used in CSPecon .round(2)
# EcombiEW = integrate.trapz(combiEW).round(2) # used in CSPecon

# capital_combiEW = (A_helio_optEW+area_combiEW)*csp_area_costs + PcombiEW*power_block_cost
# revenue_combiEW = cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW)
# cash_flow_combiEW = [-capital_combiEW] + [revenue_combiEW for i in range(30)]
# dpb_combiEW = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiEW)
# npv_combiEW = npf.npv(csp_discount_rate, [-capital_combiEW] 
# + [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])
# irr_combiEW = npf.irr([-capital_csp] 
# + [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])

# combi_finance = pd.DataFrame((dpb_combiNS,dpb_combiEW,npv_combiNS,npv_combiEW,irr_combiNS,irr_combiEW)).round(2)

# # %%

# %%
