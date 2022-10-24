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

@pytest.fixture
def sotr():
    """SolarTrough ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)
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
    pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
    Ib = pvgis.loc[:,'Gb(n)']
    return Ib

def test_perform_calc_gen_EW(sotr, Ib):
    OTrou = sotr.perform_calc(alignment='EW', Ib=Ib)
    vals = OTrou.data.values
    assert vals.sum() == pytest.approx(52583.642890174386) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(130955.74198433738 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(20.21308794388492,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(14.87011024316919,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(64.98750475959784,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.09238010496977692,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(52591.12, 1e-1) , 'Energy_MWh is incorrect'


def test_perform_calc_gen_NS(sotr, Ib):
    OTrou = sotr.perform_calc(alignment='NS', Ib=Ib)
    vals = OTrou.data.values
    assert vals.sum() == pytest.approx(74596.35065305409) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(151041.27571630833 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(21.388075760449134,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(15.253613788170941,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(62.29854949627815,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.13670332645996,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(74603.83,1e-1) , 'Energy_MWh is incorrect'


def test_perform_calcs_EW(sotr, Ib):
    OTrou = sotr.perform_calcs_EW(Ib)
    vals = OTrou.data.values
    assert vals.sum() == pytest.approx(52583.642890174386) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(130955.74198433738 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(20.21308794388492,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(14.87011024316919,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(64.98750475959784,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.09238010496977692,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(52591.12, 1e-1) , 'Energy_MWh is incorrect'


def test_perform_calcs_NS(sotr, Ib):
    OTrou = sotr.perform_calcs_NS(Ib)
    vals = OTrou.data.values
    assert vals.sum() == pytest.approx(74596.35065305409) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(151041.27571630833 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(21.388075760449134,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(15.253613788170941,1e-6) , 'Std is incorrect'
    
    assert OTrou.PowerMax_MW == pytest.approx(62.29854949627815,1e-6) , 'Max Power is incorrect'
    assert OTrou.CF == pytest.approx(0.13670332645996,1e-6) , 'CF is incorrect'
    assert OTrou.Energy_MWh == pytest.approx(74603.83,1e-1) , 'Energy_MWh is incorrect'





def test_find_units_for_maxMW(sotr, Ib):
    """ Test find_area_for_maxMW functionality
    
    Args:
        st2 (_type_): _description_
        Ib (_type_): _description_
    """    
    resEW = sotr.find_units_for_max_MW(target_MW= 50, alignment='EW', Ib= Ib,Tr=318)
    assert resEW == pytest.approx(1384.88, abs=0.05), "there is a problem with the optimisation"

    resNS = sotr.find_units_for_max_MW(target_MW= 50, alignment='NS', Ib= Ib,Tr=318)
    assert resNS == pytest.approx(1444.65, abs=0.05), "there is a problem with the optimisation"

