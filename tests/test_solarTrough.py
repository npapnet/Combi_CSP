# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""

import pathlib
import pytest 
import numpy as np
import pandas as pd 


from CombiCSP import SolarTroughCalcs, SolarSystemLocation

from datetime import datetime
def custom_date_parser(date_str):
    return datetime.strptime(date_str, "%Y%m%d:%H%M")  # Adjust format as needed

@pytest.fixture
def sotr():
    """SolarTrough ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24,  dt_gmt_hr=+2, alt=0)
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
def Ib():
    """Solar Trough ExampleData
    """
    # pytest is configured in vscode at the root
    FNAME = pathlib.Path('tests/example_data/tmy_35.015_25.755_2005_2020.csv')
    pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, 
            parse_dates=['time(UTC)'], 
            date_format = custom_date_parser, 
            engine='python') 
    Ib = pvgis.loc[:,'Gb(n)']
    return Ib

def test_perform_calc_gen_EW(sotr, Ib):
    OTrou = sotr.perform_calc(alignment='EW', Ib=Ib)
    vals = OTrou.power_data_as_array()
    assert vals.sum() == pytest.approx(52583.642890174386) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(130955.74198433738 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(20.21308794388492,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(14.87011024316919,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(64.98750475959784,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.09238010496977692,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(52591.12, 1e-1) , 'Energy_MWh is incorrect'


def test_perform_calc_gen_NS(sotr, Ib):
    OTrou = sotr.perform_calc(alignment='NS', Ib=Ib)
    vals = OTrou.power_data_as_array()
    assert vals.sum() == pytest.approx(74596.35065305409) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(151041.27571630833 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(21.388075760449134,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(15.253613788170941,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(62.29854949627815,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.13670332645996,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(74603.83,1e-1) , 'Energy_MWh is incorrect'


def test_perform_calcs_EW(sotr, Ib):
    OTrou = sotr.perform_calc(Ib=Ib, alignment='EW')
    vals = OTrou.power_data_as_array()
    assert vals.sum() == pytest.approx(52583.642890174386) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(130955.74198433738 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(20.21308794388492,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(14.87011024316919,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(64.98750475959784,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.09238010496977692,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(52591.12, 1e-1) , 'Energy_MWh is incorrect'


def test_perform_calcs_NS(sotr, Ib):
    OTrou = sotr.perform_calc(Ib=Ib, alignment='NS')
    vals = OTrou.power_data_as_array()
    assert vals.sum() == pytest.approx(74596.35065305409) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(151041.27571630833 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(21.388075760449134,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(15.253613788170941,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(62.29854949627815,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.13670332645996,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(74603.83,1e-1) , 'Energy_MWh is incorrect'



def test_incident_energy_NS(sotr, Ib):
    vals = sotr.incident_energy_on_system(Ib=Ib, alignment='NS')
    expected = [1981411.5667823767, 1981411.5667823767, 302.48040170206184  ,302.48040170206184 ]
    assert vals.sum() == pytest.approx(expected[0] ) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(expected[1],1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(expected[2],1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(expected[3],1e-6) , 'Std is incorrect'


def test_incident_energy_EW(sotr, Ib):
    vals = sotr.incident_energy_on_system(Ib=Ib, alignment='EW')
    expected = [1670115.0899976315 , 1670115.0899976315 ,285.8631617628524  ,285.8631617628524 ]
    assert vals.sum() == pytest.approx(expected[0] ) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(expected[1],1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(expected[2],1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(expected[3],1e-6) , 'Std is incorrect'



def test_find_units_for_maxMW(sotr, Ib):
    """ Test find_area_for_maxMW functionality
    
    Args:
        st2 (_type_): _description_
        Ib (_type_): _description_
    """
    expected = 1385 #1384.88
    resEW = sotr.find_units_for_max_MW(target_MW= 50, alignment='EW', Ib= Ib,Tr=318)
    assert resEW == pytest.approx(expected , abs=0.05), "there is a problem with the optimisation"

    expected = 1445 #1444.65
    resNS = sotr.find_units_for_max_MW(target_MW= 50, alignment='NS', Ib= Ib,Tr=318)
    assert resNS == pytest.approx(expected, abs=0.05), "there is a problem with the optimisation"


def test_calculate_nominal_power(sotr):
    """ Test calculate_nominal_power functionality
    
    Args:
        sotr (_type_): _description_
    """
    expected = 63.227
    res = sotr.calculate_nominal_power_MW(Tr=318, T_amb=15, nG=0.97)
    assert res == pytest.approx(expected, abs=1e-1), "Nominal power is incorrect"

    
def test_find_no_of_units_for_nominalPower (sotr):
    """ Test calculate_nominal_power functionality
    
    Args:
        sotr (_type_): _description_
    """
    expected = 1424
    res = sotr.find_no_units_for_nominal_power_MW(target_MW=50, Tr=318)
    assert res == pytest.approx(expected, abs=1e-1), "Nominal power is incorrect"