# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""

import pytest 
import numpy as np
from CombiCSP.solar_system_location import SolarSystemLocation
import CombiCSP.solar_system_location as ccspSSL

@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])
@pytest.fixture
def Crete():
    """Example site
    """    
    return SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt_hr=+2, alt=0)

@pytest.fixture
def R():
    """R for transmittance functions 
    units are unknown. For verifications see : 
    H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134."
    """
    return 1

    
     
class Test_AM:
    """testing the air mass functions
    """    
    def test_AM(self, hoy_ex1):
        
        # assert sgh.AM(hoy_ex1)[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        # assert sgh.air_mass(hoy_ex1, method='wiki')[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        # assert sgh.air_mass(hoy_ex1 )[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        # pytest.approx(sgh.AM(hoy_ex1, method='pveducation')[0], -13.37218176, abs=1e-3)
        pass

    def test_AM2(self,hoy_ex1):
        # assert sgh.AM2(hoy_ex1)[0] == pytest.approx( -1.1149391977788454, abs=1e-3)
        # assert sgh.air_mass(hoy_ex1, method='Kasten')[0] == pytest.approx( -1.1149391977788454, abs=1e-3)
        pass
        
    def test_AM3(self,hoy_ex1 ):
        # assert sgh.AM3(hoy_ex1)[0] == pytest.approx( -1.1184591949819083, abs=1e-3)
        # assert sgh.air_mass(hoy_ex1, method='Kasten-Young')[0] == pytest.approx( -1.1184591949819083, abs=1e-3)
        pass

    def test_AM4(self,hoy_ex1 ):
        # assert sgh.AM4(hoy_ex1)[0] == pytest.approx( 1374.4966167568873, abs=1e-3)
        # assert sgh.air_mass(hoy_ex1, method='Schoenberg')[0] == pytest.approx( 1374.4966167568873, abs=1e-3)
        pass
        


class Test_EOT_zen_ele:
    """testing the air mass functions
    """    
    def test_EOT(self, hoy_ex1, Crete):
        expected = -2.9041689600000002 
        assert ccspSSL.EoT(hoy_ex1)[0] == pytest.approx(  expected , abs=1e-3)

         
    def test_tsol(self,hoy_ex1, Crete):
        expected = 0.551597184 
        # assert sgh.tsol(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.tsol(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        
    def test_ele(self, hoy_ex1, Crete ):
        expected = -1.3257002183993918 
        # assert sgh.ele(hoy_ex1)[0] == pytest.approx(expected, abs=1e-3)
        assert Crete.ele_rad(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        
    def test_z(self,hoy_ex1, Crete ):
        expected =  2.8964965451942883 
        # assert sgh.z(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.z_rad(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)

    def test_azim(self,hoy_ex1, Crete):
        expected =  -0.577725020643127
        # assert sgh.azim(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.azim_rad(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
    

# def test_solar_hour_angle_W_deg(Crete):
#     Crete_offset = -7.8807413669699855
#     assert Crete.W(12) == pytest.approx(0+Crete_offset, abs=1e-3)
#     # assert Crete.W(0) == pytest.approx(180-Crete_offset, abs=1e-3)
#     assert Crete.W(24) == pytest.approx(180-Crete_offset, abs=1e-3)
#     assert Crete.W(6) == pytest.approx(-90+Crete_offset, abs=1e-3)
#     assert Crete.W(11) == pytest.approx(-15+Crete_offset, abs=1e-3)

