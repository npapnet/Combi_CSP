# -*- coding: utf-8 -*- 
"""
    @Author: G. Arnaoutakis, N. Papadakis
    @Date: 2022/  
"""

#TODO this module contains unused functions (not  tested)
#TODO    theta_transversal  function exists twice
#%%
'''Concentrating Solar Power plants'''
import ssl
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from iapws import IAPWS97
from CombiCSP import CtoK, HOYS_DEFAULT, SolarSystemLocation

#%%

#%%
# generator data

# SM = 2.5 # Solar Multiplier


#%% ===========================================================================
def theta_transversal(ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT)->np.array : 
    """Parabolic Trough theta  transversal incidence angle

    #TODO This function has the same name with another one in the same module
    
    #TODO  not tested
    

    Buscemi, A.; Panno, D.; Ciulla, G.; Beccali, M.; Lo Brano, V. 
    Concrete Thermal Energy Storage for Linear Fresnel Collectors: 
    Exploiting the South Mediterranean’s Solar Potential for Agri-Food Processes. 
    Energy Conversion and Management 2018, 166, 719–734, doi:10.1016/j.enconman.2018.04.075.

    Args:
        hoy (np.array): hour of year 

    Returns:
        float: theta  transversal incidence angle
    """    

    return np.arctan(np.sin(np.radians(ssloc.azim(hoy))) * np.tan(np.radians(ssloc.z(hoy))))


def theta_i(ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT)->float: 
    """Parabolic Trough longitudinal incidence angle

    Buscemi, A.; Panno, D.; Ciulla, G.; Beccali, M.; Lo Brano, V. 
    Concrete Thermal Energy Storage for Linear Fresnel Collectors: 
    Exploiting the South Mediterranean’s Solar Potential for Agri-Food Processes. 
    Energy Conversion and Management 2018, 166, 719–734, doi:10.1016/j.enconman.2018.04.075.
    
    #TODO  not tested

    Args:
        hoy (np.array): hour of year 

    Returns:
        float: not tested
    """    
    return np.arctan(np.cos(np.radians(ssloc.azim(hoy))) * np.tan(np.radians(ssloc.z(hoy)))* np.cos(theta_transversal(ssloc=ssloc)))


# not tested
def thetai_transversal(ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT):
    """_summary_

    #TODO This function has the same name with another one in the same module
    
    Morin, G.; Dersch, J.; Platzer, W.; Eck, M.; Häberle, A. 
    Comparison of Linear Fresnel and Parabolic Trough Collector Power Plants. 
    Solar Energy 2012, 86, 1–12, doi:10.1016/j.solener.2011.06.020.

    
    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

    Returns:
        _type_: _description_
    """    
    return np.arctan(abs(np.sin(ssloc.azim(hoy)))/np.tan(ssloc.ele(hoy)))

def thetai_longtitudinal(ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT): 
    return np.arcsin(np.cos(ssloc.azim(hoy))*np.cos(ssloc.ele(hoy)))


def thetai(ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT, inclination=90, azimuths=0): # incidence angle [in radians]
    #TODO This function originially resides is solar_system_location.py
    g = np.degrees(ssloc.azim(hoy)) - azimuths # if surface looks due S then azimuths=0
    return np.arccos(np.cos(ssloc.ele(hoy)) * np.sin(np.radians(inclination)) * np.cos(np.radians(g)) 
        + np.sin(ssloc.ele(hoy)) * np.cos(np.radians(inclination)))

def shade_function(Ws,Wc, ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT):
    """Shade function for Parabolic Trough Solar Power Plants, 

    Stuetzle (2002) pp.29 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
    N. Fraidenraich, C. Oliveira, A.F. Vieira da Cunha, J.M. Gordon, O.C. Vilela, 
    Analytical modeling of direct steam generation solar power plants, Solar Energy. 98 (2013) 511–522. 
    https://doi.org/10.1016/j.solener.2013.09.037.
    
    Args:
        Ws (_type_): #TODO  solar angle
        Wc (_type_): _description_
        hoy (_type_): _description_

    Returns:
        _type_: _description_
    """    
    return abs(Ws * np.cos(ssloc.z(hoy)) / (Wc * np.cos(thetai(ssloc, hoy))))

def end_loss(f,L,N, ssloc:SolarSystemLocation, hoy:np.array=HOYS_DEFAULT):
    '''Lippke, 1995 in pp.31 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
    return (1 - (f * np.tan(thetai(ssloc=ssloc, hoy=hoy)) / L)) * N

# def loss_regr(input_dict):
#     '''pp.36-42 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
#     University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
#     '''
#     equation = a0 + a1 * T + a2 * T**2 + a3 * T**3 + Ib2(alt) * (b0 + b1 * T**2)
#     # Use dict flag to get {variable: value} output, not anonymous [value]
#     #solution = solve(equation.subs(input_dict), dict=True)
#     return equation

