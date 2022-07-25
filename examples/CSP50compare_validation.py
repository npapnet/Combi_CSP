#%% this file contains code from the original code and the new OOP code, and perfoms in line comparison
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

from CombiCSP.solar_tower import solarII,IAM_tow
from CombiCSP.solar_trough import di_sst, IAM_tro, costhetai_NS, costhetai_EW, Ac, Cg_tro

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
# capital_csp = 5000000 # this is obsolete because its calculated from the characteristics of the systems


# Tower dimensions
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m]
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
#R = 5


area_list = []
cash_flow_list = []
tow_scenaria = []
for A_helio in np.arange(75000,125001,10000): # 100MW np.arange(150000,250001,10000):
    Ctow = A_helio / Ar
    tower = solarII(Ib,1,IAM_tow(hoy),A_helio,Ar)
    tow_xyz = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
    tow_mon = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
    Ptower = np.amax(tower) # used in CSPecon .round(2)
    Etower = integrate.trapz(tower).round(2) # used in CSPecon
    CF_tow = Etower / (8760 * Ptower)#.round(2)
    CF_towh = Etower / (8760 * tower) #<<<<<<<<<<<<<<<<<<<<<<<<<too noisy, test with positive data
    tow_data = np.vstack((A_helio,Ctow,Ptower,Etower,CF_tow)) # vertical stack
    # economics
    capital_csp_tow = A_helio*csp_area_costs + Ptower*power_block_cost
    revenue_csp_tow = cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow)
    cash_flow_tow = [-capital_csp_tow] + [revenue_csp_tow for i in range(30)]
    dpb_tow = cspe.discounted_payback_period(csp_discount_rate, cash_flow_tow)
    npv_csp_tow = npf.npv(csp_discount_rate, [-capital_csp_tow] + [cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow) for i in range(30)])
    irr_csp_tow = npf.irr([-capital_csp_tow] + [cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow) for i in range(30)])
    area_list.append(A_helio)
    cash_flow_list.append(cash_flow_tow)
    tow_scenaria.append((A_helio,Ctow,Ptower,Etower,CF_tow,dpb_tow,npv_csp_tow,irr_csp_tow,cash_flow_tow))
    # plt.plot(hoy, tower, label=A_helio)
# plt.xlabel('Time (hour of year)')
# plt.ylabel('Power (MW)') 
# plt.title('Tower')
# plt.legend()
# #xlim(0,87.60), ylim(0,80)
# plt.show()



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
            # capital_csp=capital_csp,
        lifetime=range(30))
    area_list3.append(tmp_res_Dic['A_helio'] )
    cash_flow_list3.append(tmp_res_Dic['cash_flow'])
    tow_scenaria3.append(tmp_res_Dic['scenaria'])



#%% Assertions 
np.testing.assert_equal(area_list, area_list3) 
np.testing.assert_equal(cash_flow_list, cash_flow_list3) 
for k in range(len(tow_scenaria3)):
    np.testing.assert_equal(tow_scenaria3[k][:8], tow_scenaria[k][:8]) 
    np.testing.assert_equal(tow_scenaria3[k][8], tow_scenaria[k][8]) 

print('tests between original code and OOP completed without problems')

#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Wr = 0.07 # tube outer diameter [m]
Wc = 5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine


## the following line should be uncommented. Otherwise the lists will grow from the Tower.
# area_list = []  
# cash_flow_list = []
trough_scenaria = []
troughew_scenaria = []
for N in np.arange(800,1301,100): # 100MW np.arange(1000,2001,100):
    area = Ac(Wc, L, N)
    # NS calcs
    trough = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
    tro_xyz = np.vstack(trough).reshape((365,24)) # reshape 8760,1 to 365,24
    Ptrough = np.amax(trough) # used in CSPecon .round(2)
    Etrough = integrate.trapz(trough).round(2) # used in CSPecon
    CF_tro = Etrough / (8760 * Ptrough)#.round(2)
    tro_data = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro)) # vertical stack
    
    # EW calcs
    troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
    datah = np.vstack((hoy, trough))
    Ptroughew = np.amax(troughew) # used in CSPecon .round(2)
    Etroughew = integrate.trapz(troughew).round(2) # used in CSPecon
    CF_troew = Etroughew / (8760 * Ptroughew)#.round(2)
    tro_dataew = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew)) # vertical stack

    #plots
    plt.plot(hoy, trough, label=N)#,xlim(100,600)
    plt.plot(hoy, troughew, label=N)#,xlim(100,600)
    # economics
    capital_csp_tro = area*csp_area_costs + Ptrough*power_block_cost
    revenue_csp_tro = cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro)
    cash_flow_tro = [-capital_csp_tro] + [revenue_csp_tro for i in range(30)]
    dpb_tro = cspe.discounted_payback_period(csp_discount_rate, cash_flow_tro)
    npv_csp_tro = npf.npv(csp_discount_rate, [-capital_csp_tro] \
        + [cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro) for i in range(30)])
    irr_csp_tro = npf.irr([-capital_csp_tro] \
        + [cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro) for i in range(30)])
    area_list.append(area)
    cash_flow_list.append(cash_flow_tro)
    trough_scenaria.append((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro,dpb_tro,npv_csp_tro,irr_csp_tro,cash_flow_tro))
    troughew_scenaria.append((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew,dpb_tro,npv_csp_tro,irr_csp_tro,cash_flow_tro))
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)'), 
plt.legend()
plt.title('Trough')
#xlim(0,87.60), ylim(0,80)
plt.show()

#%%

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

#%% Assertions 
np.testing.assert_equal(area_list[6:], ns_area_list) 
np.testing.assert_equal(cash_flow_list[6:], ns_cash_flow_list) 
for k in range(len(ns_trough_scenaria)):
    np.testing.assert_equal(ns_trough_scenaria[k][:8], trough_scenaria[k][:8]) 
    np.testing.assert_equal(ns_trough_scenaria[k][8], trough_scenaria[k][8]) 

print('tests between original code and OOP completed without problems')
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

DPB = [] #discounted payback period
for y in cash_flow_list:
    DPB.append(cspe.discounted_payback_period(csp_discount_rate, y).round(2))
dpb_table = pd.DataFrame(DPB, index=area_list, columns = ['DPB (years)'])
dpb_table

# T+NS optimum
A_helio_optNS = 125000
N_opt_NS = 800  # TODO find out where does this come from

#=== NS optimisation
tower_opt = solarII(Ib,1,IAM_tow(hoy),A_helio_optNS,Ar)
trough_opt = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_NS)
combiNS = tower_opt + trough_opt
combiNS_xyz = np.vstack(combiNS).reshape((365,24)) # reshape 8760,1 to 365,24
plt.title('Tower + Trough N-S')
cspm.heatmap2d(combiNS_xyz.T)

area_combiNS = Ac(Wc, L, N_opt_NS)
PcombiNS = np.amax(combiNS) # used in CSPecon .round(2)
EcombiNS = integrate.trapz(combiNS).round(2) # used in CSPecon

capital_combiNS = (A_helio_optNS+area_combiNS)*csp_area_costs + PcombiNS*power_block_cost
revenue_combiNS = cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS)
cash_flow_combiNS = [-capital_combiNS] + [revenue_combiNS for i in range(30)]
dpb_combiNS = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiNS)
npv_combiNS = npf.npv(csp_discount_rate, [-capital_combiNS] \
    + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])
irr_combiNS = npf.irr([-capital_combiNS] \
    + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])


#%%
stc_opt =  SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, 
        Ar = 99.3 , A_helio = A_helio_optNS,
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

#%% Assertions 
np.testing.assert_equal(combiNS, combiNS_np) 
np.testing.assert_equal(combiNS_xyz, combiNS_xyz_np) 
assert area_combiNS == oTr.A_helio
assert PcombiNS == (oTr.data + oTow.data).max()
assert EcombiNS == (oTr.Energy_MWh + oTow.Energy_MWh).round(2)

assert (A_helio_optNS+area_combiNS) == tmp_res_Dic['A_helio']
np.testing.assert_equal( cash_flow_combiNS,  tmp_res_Dic['cash_flow'])
assert (A_helio_optNS+area_combiNS) ==  tmp_res_Dic['scenaria'][0]
assert   tmp_res_Dic['scenaria'][1] is None
assert   PcombiNS == tmp_res_Dic['scenaria'][2] , 'Max Power '
assert   EcombiNS == tmp_res_Dic['scenaria'][3].round(2) , 'Total Energy Yield '
assert   tmp_res_Dic['scenaria'][4] is None , 'Capacity Factor '
assert   tmp_res_Dic['scenaria'][5] == dpb_combiNS , 'discounted payback '
assert   tmp_res_Dic['scenaria'][6] == npv_combiNS , 'net present value'
assert   tmp_res_Dic['scenaria'][7] == irr_combiNS , 'internal rate of return'
np.testing.assert_equal( cash_flow_combiNS,  tmp_res_Dic['scenaria'][8])
print('tests between original code and OOP completed without problems')
#%% ========================================================================
# T+EW optimum
A_helio_optEW = 125000
N_opt_EW = 800
tower = solarII(Ib,1,IAM_tow(hoy),A_helio_optEW,Ar)
troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_EW)
combiEW = tower + troughew
combiEW_xyz = np.vstack(combiEW).reshape((365,24)) # reshape 8760,1 to 365,24
plt.title('Tower + Trough E-W')
cspm.heatmap2d(combiEW_xyz.T)

area_combiEW = Ac(Wc, L, N_opt_EW)
PcombiEW = np.amax(combiEW) # used in CSPecon .round(2)
EcombiEW = integrate.trapz(combiEW).round(2) # used in CSPecon

capital_combiEW = (A_helio_optEW+area_combiEW)*csp_area_costs + PcombiEW*power_block_cost
revenue_combiEW = cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW)
cash_flow_combiEW = [-capital_combiEW] + [revenue_combiEW for i in range(30)]
dpb_combiEW = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiEW)
npv_combiEW = npf.npv(csp_discount_rate, [-capital_combiEW] 
+ [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])
irr_combiEW = npf.irr([-capital_combiEW] 
+ [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])

combi_finance = pd.DataFrame((dpb_combiNS,dpb_combiEW,npv_combiNS,npv_combiEW,irr_combiNS,irr_combiEW)).round(2)

# %%

# %%
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

#%% Assertions 
np.testing.assert_equal(combiEW, combiEW_np) 
np.testing.assert_equal(combiEW_xyz, combiEW_xyz_np) 
assert area_combiNS == oTr.A_helio
assert PcombiNS == (oTr.data + oTow.data).max()
assert EcombiNS == (oTr.Energy_MWh + oTow.Energy_MWh).round(2)

assert (A_helio_optNS+area_combiNS) == tmp_res_Dic['A_helio']
np.testing.assert_equal( cash_flow_combiEW,  tmp_res_Dic['cash_flow'])
assert (A_helio_optNS+area_combiNS) ==  tmp_res_Dic['scenaria'][0]
assert   tmp_res_Dic['scenaria'][1] is None
assert   PcombiEW == tmp_res_Dic['scenaria'][2] , 'Max Power '
assert   EcombiEW == tmp_res_Dic['scenaria'][3].round(2) , 'Total Energy Yield '
assert   tmp_res_Dic['scenaria'][4] is None , 'Capacity Factor '
assert   tmp_res_Dic['scenaria'][5] == dpb_combiEW , 'discounted payback '
assert   tmp_res_Dic['scenaria'][6] == npv_combiEW , 'net present value'
assert   tmp_res_Dic['scenaria'][7] == irr_combiEW , 'internal rate of return'
np.testing.assert_equal( cash_flow_combiEW,  tmp_res_Dic['scenaria'][8])
print('tests between original code and OOP completed without problems')
# %%
