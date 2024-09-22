# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""
#NOTE these are fictional tests used in order to get a base value. They should be properly validated.
import pytest 
import numpy as np

import CombiCSP.CSP as cspC
from CombiCSP import SolarSystemLocation

@pytest.fixture
def Crete():
    """Example site
    """    
    return SolarSystemLocation(lat=35, lon=24,  dt_gmt_hr=+2, alt=0)


@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])

def test_theta_transversal(hoy_ex1, Crete):
    expected = -0.0005101662427872909
    assert cspC.theta_transversal(ssloc=Crete, hoy =hoy_ex1)[0] == pytest.approx(expected, abs=1e-6)
    
    actual = cspC.theta_transversal(ssloc=Crete)
    assert  actual.sum() == pytest.approx(   8.72523757357484e-05, abs = 1e-6 )
    assert  actual.std() == pytest.approx(  0.0005363889338115103, abs = 1e-6 )


def test_theta_i(hoy_ex1, Crete):
    expected = 0.05055082940933891
    assert cspC.theta_i(ssloc=Crete, hoy=hoy_ex1)[0] == pytest.approx(expected, abs=1e-6)
    
    actual = cspC.theta_i(ssloc=Crete, hoy=hoy_ex1)
    assert  actual.sum() == pytest.approx(  442.8252595598201, abs = 1e-6 )
    assert  actual.std() == pytest.approx(  7.38237756908075e-09, abs = 1e-6 )


def test_thetai_longtitudinal(hoy_ex1, Crete):
    expected = 0.204695771081944
    assert cspC.thetai_longtitudinal(ssloc=Crete, hoy=hoy_ex1)[0] == pytest.approx(expected, abs=1e-6)
    
    actual = cspC.thetai_longtitudinal(ssloc=Crete, hoy=hoy_ex1)
    assert  actual.sum() == pytest.approx(  0.204695771081944, abs = 1e-6 )
    assert  actual.std() == pytest.approx(  0.0, abs = 1e-6 )


def test_shade_function(hoy_ex1, Crete):
    Ws = 1
    Wc = 2
    expected = 2.3862778289930486
    assert cspC.shade_function(Ws,Wc,ssloc=Crete, hoy=hoy_ex1)[0] == pytest.approx(expected, abs=1e-6)
    
    actual = cspC.shade_function(Ws,Wc,ssloc=Crete)
    assert  actual.sum() == pytest.approx(  30012.324364512737, abs = 1e-6 )
    assert  actual.std() == pytest.approx(  34.882632160821885, abs = 1e-6 )

def test_end_loss(hoy_ex1, Crete):
    f = 0.88
    L = 25
    N = 1800
    Ws = 1
    Wc = 2
    expected = 1494.8027687770161
    assert cspC.end_loss(f,L,N,ssloc=Crete, hoy=hoy_ex1)[0] == pytest.approx(expected, abs=1e-6)
    
    actual = cspC.end_loss(f,L,N,ssloc=Crete)
    assert  actual.sum() == pytest.approx(  8211925.231453498, abs = 1e-6 )
    assert  actual.std() == pytest.approx(  9110.752744711397, abs = 1e-6 )
