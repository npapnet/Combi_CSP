import numpy as np
import pandas as pd
from CombiCSP.financial.lcoe import LCOECalculator

import pytest

@pytest.fixture
def df_annual_data():
    year  =  range(0, 31)
    energy_MWh = np.ones(31) * (1600)
    costs_per_annum = np.ones(31) * 10000
    costs_per_annum[0] = 1e6

    df = pd.DataFrame({'year': year,
        'energy_MWh':energy_MWh, 
        'costs_per_annum': costs_per_annum
        })    
    return df
    

def test_lcoe_calculator( df_annual_data):
    
    year = df_annual_data['year'].values
    discount_rate = 0.04
    
    discounted_costs = df_annual_data['costs_per_annum'] / (1 + discount_rate) ** year
    
    disc_factor = 1 / (1 + discount_rate) ** year
    
    np.testing.assert_array_almost_equal(discounted_costs[1:], 10000*disc_factor[1:])
        
    lcoecalc = LCOECalculator() 
    expected = 40.07
    actual = lcoecalc.lcoe_Eu_MWh(df=df_annual_data, discount_rate=discount_rate)
    assert expected == pytest.approx(actual, abs=1e-2)
    