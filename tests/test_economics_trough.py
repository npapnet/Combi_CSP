import pathlib
import numpy as np
import pandas as pd
import pytest
from CombiCSP import EconomicEnvironment, SolarSystemLocation,  SolarTroughCalcs, HOYS_DEFAULT

from datetime import datetime
def custom_date_parser(date_str):
    return datetime.strptime(date_str, "%Y%m%d:%H%M")  # Adjust format as needed

@pytest.fixture
def stro()->SolarTroughCalcs:
    """SolarTrough ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24, dt_gmt_hr=+2, alt=0)
    return SolarTroughCalcs(
        foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
        ,N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
        ,L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  
        ,Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
        ,Wr = 0.07 # tube outer diameter [m]
        ,Wc = 5.76,
        slobj= slobj
        )


@pytest.fixture
def ee():
    """Economic environment 
    """    
    return EconomicEnvironment(   
            oil_price=60, 
            Eoil=11.766,
            currency_units='USD')
    
@pytest.fixture
def Ib():
    """SolarTower ExampleData
    """
    # pytest is configured in vscode at the root
    FNAME = pathlib.Path('tests/example_data/tmy_35.015_25.755_2005_2020.csv')
    pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, 
            parse_dates=['time(UTC)'], 
            date_format = custom_date_parser,
            engine='python') 
    Ib = pvgis.loc[:,'Gb(n)']
    return Ib
    


def test_trough_NS_fin(stro, ee, Ib):
    o_tmp = stro.perform_calc(alignment='NS', Ib= Ib, Tr=318, hoy=HOYS_DEFAULT)
    tmp_res_Dic = stro.financial_assessment(oTr=o_tmp, 
        ee = ee,
        csp_area_costs = 235, 
        csp_energy_price = 248, 
        csp_discount_rate = 0.09, 
        power_block_cost = 910000.0,
        lifetime=30
        )
    
    scenario_dic  = tmp_res_Dic.get('scenario_params', None)
    assert scenario_dic.get('area_m2') == 259200.0 
    assert scenario_dic.get('system_type') == 'trough_NS'


    cash_flow  = tmp_res_Dic['cash_flow_df']['cash_flow'].values
    np.testing.assert_almost_equal(
            cash_flow[:2],
            [-117603680.04161312,   13795837.74])


    syst_res = tmp_res_Dic['scenario_params']
    fin_res = tmp_res_Dic['scenario_financial']

    expecting = (259200.0, # area_m2
        82.28571428571428, # Ctow
        62.29854949627815, # PowerMax_MW
        74603.83,          # Energy_MWh
        0.13670332645996,  # CF
        16.91753934467728, # discounted_payback_period
        24129984.13280969, # npv
        0.11252063804209933) # irr
    
    assert syst_res.get('Cg') == expecting[1], 'Cg - concentration ratio'

    assert fin_res.get('PowerMax_MW') == pytest.approx(expecting[2], abs=1e-3), 'PowerMax_MW'
    assert fin_res.get('Energy_MWh') == pytest.approx(expecting[3], abs=1e-3), 'Energy_MWh'
    assert fin_res.get('CF') == pytest.approx(expecting[4], abs=1e-3), 'CF'
    assert fin_res.get('discounted_payback_period') == pytest.approx(expecting[5], abs=1e-3), 'discounted_payback_period'
    assert fin_res.get('npv') == pytest.approx(expecting[6], rel=1e-5), 'NPV'
    assert fin_res.get('irr') == pytest.approx(expecting[7], abs=1e-3), 'IRR'


def test_trough_EW_fin(stro, ee, Ib):
    o_tmp = stro.perform_calc(alignment='EW', Ib= Ib, Tr=318, hoy=HOYS_DEFAULT)
    tmp_res_Dic = stro.financial_assessment(oTr=o_tmp, 
        ee = ee,
        csp_area_costs = 235, 
        csp_energy_price = 248, 
        csp_discount_rate = 0.09, 
        power_block_cost = 910000.0,
        lifetime=30
        )
    
    


    cash_flow  = tmp_res_Dic['cash_flow_df']['cash_flow'].values
    np.testing.assert_almost_equal(
            cash_flow[:2],
            [-120050629.33123404,   8238807.69])

    syst_res = tmp_res_Dic['scenario_params']
    assert syst_res.get('system_type') == 'trough_EW'
    assert syst_res.get('area_m2') == 259200.0 
    assert syst_res.get('Ac_m2') == 259200.0 
    assert syst_res.get('Ar_m2') == pytest.approx(3150, abs=1e-3) 
    assert syst_res.get('Cg') == pytest.approx(82.28571428571428, abs=1e-3), 'Cg - concentration ratio'


    fin_res = tmp_res_Dic['scenario_financial']

    expecting = (259200.0, # area_m2
        82.28571428571428, # Ctow
        64.98750475959784, # PowerMax_MW
        52591.12,          # Energy_MWh
        0.09238010496977692,  # CF
        -1, # discounted_payback_period
        -35407969.39718697, # npv
        0.054766137547524796) # irr
    


    assert fin_res.get('PowerMax_MW') == pytest.approx(expecting[2], abs=1e-3), 'PowerMax_MW'
    assert fin_res.get('Energy_MWh') == pytest.approx(expecting[3], abs=1e-3), 'Energy_MWh'
    assert fin_res.get('CF') == pytest.approx(expecting[4], abs=1e-3), 'CF'
    assert fin_res.get('discounted_payback_period') == pytest.approx(expecting[5], abs=1e-3), 'discounted_payback_period'
    assert fin_res.get('npv') == pytest.approx(expecting[6], rel=1e-5), 'NPV'
    assert fin_res.get('irr') == pytest.approx(expecting[7], abs=1e-3), 'IRR'