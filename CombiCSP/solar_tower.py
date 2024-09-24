# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
    @Credit: original functions from G. Arnaoutakis
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import numpy_financial as npf

from CombiCSP import OutputContainer, CtoK, HOYS_DEFAULT, SolarSystemLocation
from CombiCSP.solar_system_location import SolarSystemLocation

from .financial.economics import cashflow, discounted_payback_period, EconomicEnvironment

  
STEFAN_BOLTZMANN_CONSTANT = 5.67 * 1e-8 # Stefan – Boltzman constant [W/m2K4]

class SolarTowerCalcs():
    _system_type = 'tower'
    results_container : OutputContainer = None
    optical_efficiency:float= 100 # heliostat optical effiency [%] 65% pp.24 in Pacheco
    
    nR:float = 0.412 # SAM default
    reflectivity:float  = 1 # reflectivity [%] 1 if IAM is IAM_tow(hoy)
    epsilon:float  = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
    alpha:float  = 1 # absorptivity of the receiver
        
    def __init__(self, 
        alt_:float = 200*10e-3 
        , Ht_km:float = 0.1
        , Ar_m2:float = 99.3 
        , A_helio_m2:float = 225000
        , slobj:SolarSystemLocation = None
        ):
        self.Ar_m2 = Ar_m2# receiver area [m2] pp.44 in Pacheco
        self.alt_km = alt_ #Height above sea level (probably [km] uncertain TODO check with original code)
        self.Ht_km = Ht_km # Tower height [km]
        #A_helio = 71140 + 10260 # total heliostat area [m2] pp.22 in Pacheco
        self.A_helio_m2 = A_helio_m2 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project
        self.Ctow = self.A_helio_m2 / self.Ar_m2
        if slobj is None:
            raise ValueError('System location not found')
        else:
            self._sl = slobj

    def params_as_dict(self)->dict:
        """Returns the parameters of the solar tower as a dictionary	
        """
        return {
            'alt_km':self.alt_km, 
            'Ht_km':self.Ht_km, 
            'Ar_m2':self.Ar_m2, 
            'A_helio_m2':self.A_helio_m2,
            'Ctow':self.Ctow
            , 'optical_efficiency':self.optical_efficiency
            , 'reflectivity':self.reflectivity
            , 'nR':self.nR
            , 'emissivity':self.epsilon
            , 'absorptivity':self.alpha
            }
    
    def perform_calc(self, Ib:pd.Series, transmittance=1, nG:float=0.97, hoy:np.array=HOYS_DEFAULT)->OutputContainer:
        """Performs solar tower calculations

        Args:
            Ib (_type_): _description_
            transmittance (int, optional): _description_. Defaults to 1.
            nG (float): Generator efficiency (TODO Crosscheck??)
            hoy (_type_, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            _type_: _description_
        """
        # data = solarII(Ib=Ib,Trans=transmittance, IAM=IAM_tow(hoy)
        
        hourly_energy_yield = self.solarII(Ib=Ib,Trans=transmittance, nG=0.97, hoy=hoy)
        df = pd.DataFrame({'HOY':hoy,'Ib_n':Ib, 'Power_MW':hourly_energy_yield})
        tower_params = {'system_type':self._system_type} 
        tower_params.update(self.params_as_dict())
        
        self.results_container = OutputContainer(
                power_df = df,
                scenario_params = tower_params,
                system_type= self._system_type
            )
        return self.results_container
    
    def solarII(self, Ib:pd.Series,Trans:float, nG:float = 0.97, hoy:np.array=HOYS_DEFAULT
            , Trec_Celcius:float = 565 # the working fluid temperature in the receiver [oC] 565 oC 838K
            , Ta_Celcius:float = 15 # ambient temperature close to the receiver [oC] 15 oC 288K
            , Tin_Celcius:float = 290 # working fluid inlet temperature to the receiver [oC] 290 oC 563K
            )->pd.Series:
        """Calculates the power of the solar tower with heliostat
    
        R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
        Energy Conversion and Management. 60 (2012) 226–232.
        J.E. Pacheco, R.W. Bradshaw, D.B. Dawson, W.D. la Rosa, R. Gilbert, S.H. Goods, P. Jacobs, M.J. Hale, S.A. Jones, G.J. Kolb, M.R. Prairie, H.E. Reilly, S.K. Showalter, L.L. VANT-HULL, 
        Final Test and Evaluation Results from the Solar Two Project, n.d. https://core.ac.uk/reader/193342950 (accessed September 8, 2020).
    
        Args:
            Ib (pd.Series): direct irradiance
            Trans (float): transmissivity 
            IAM (np.array): incidence angle modifier
            A_helio (float): heliostat area in m^2
            Ar (float): receiver area in m^2
            nG (float):  efficiency of generator Mosleh19 (Defaults to 0.97)

        Returns:
            pd.Series: power in MW (since the output is hourly, this can be considered as energy in MWh)
        """
        IAM  = self.IAM_tow(hoy=hoy)
        Heliostat_area = self.A_helio_m2
        ReceiverArea = self.Ar_m2
        
        # Incident energy on the receiver
        Qin = Ib * self.reflectivity * Trans * IAM * Heliostat_area * self.optical_efficiency/100 # Eq. 17  in McGovern12
        
        # Thermal energy radiated by the receiver
        Qrad = self.epsilon * STEFAN_BOLTZMANN_CONSTANT * ReceiverArea * (CtoK(Trec_Celcius)**4-CtoK(Ta_Celcius)**4)

        # Convection losses
        hconv = CtoK(Trec_Celcius)/60 + 5/3 # [W/m2K] convection coefficient Eq. 20 in McGovern12
        '''D.L. Siebers, J.S. Kraabel, Estimating convective energy losses from solar central 
        receivers, Sandia National Lab. (SNL-CA), Livermore, CA (United States), 1984. 
        https://doi.org/10.2172/6906848.'''
        Qconv = hconv * ReceiverArea * (CtoK(Trec_Celcius) - CtoK(Ta_Celcius))
        
        # Energy balance
        Qnet = self.alpha * Qin - Qrad - Qconv # Eq. 8 in McGovern12
        
        if Qnet.all() <= 0: #<<<<<<<<<<<<<<<<<<< check with Qin
            P = 0
        else:
            P = Qnet * self.nR * nG
        return P/1e6 # convert W to MW

    # Incidence angle methods for towers

    def IAM_tow(self, hoy:np.array=HOYS_DEFAULT)->np.array : 
        """Incidence angle modifier of Tower (azimuth)

        for explanation see: http://www.solarpanelsplus.com/solar-tracking/
        
        # polynomial fit, see file IAM.py for data
        
        Args:
            hoy (np.array): hour of year

        Returns:
            np.array : Incidence angle modifier of Tower in rad
        """    
        zenith_deg = np.degrees(self._sl.z_rad(hoy))
        return 1.66741484e-1 + 1.41517577e-2 * zenith_deg - 9.51787164e-5 * zenith_deg**2
    
    def incident_energy_on_system(self,  Ib:pd.Series, hoy:np.array = HOYS_DEFAULT)->pd.Series:
        """Calculates the incident energy on the system

        Args:
            Ib (pd.Series): direct irradiance
            hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.
            
        Returns:
            pd.Series: incident energy on the system
        """
        return Ib*self.IAM_tow(hoy)

    def mutate(self,   alt = None , Ht = None
        , Ar = None 
        , A_helio = None
        , slobj:SolarSystemLocation = None):
        """Function that produces a new Solar tower object with different parameters

        If a parameter is not set then it does not change.

        Args:
            alt (_type_, optional): _description_. Defaults to None.
            Ht (_type_, optional): _description_. Defaults to None.
            Ar (_type_, optional): _description_. Defaults to None.
            A_helio (_type_, optional): _description_. Defaults to None.
            slobj (SolarSystemLocation, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """        
        alt = self.alt_km if alt is None else alt
        Ht = self.Ht_km if Ht is None else Ht
        Ar = self.Ar_m2 if Ar is None else Ar
        A_helio = self.A_helio_m2 if A_helio is None else A_helio
        slobj = self._sl if slobj is None else slobj
        return SolarTowerCalcs(alt_ = alt, Ht_km=Ht, Ar_m2=Ar, A_helio_m2=A_helio, slobj=slobj)

    def find_area_for_max_MW(self, target_MW :float, 
            Ib, 
            transmittance=1, nG=0.97, 
            hoy=HOYS_DEFAULT)->float:
        """finds the required area for the parameters of the solar tower

        Args:
            target_MW(float): max MW targed
            Ib (_type_): _description_
            transmittance (int, optional): _description_. Defaults to 1.
            nG (float, optional): _description_. Defaults to 0.97.
            hoy (_type_, optional): _description_. Defaults to HOYS_DEFAULT.
        
        returns:
            (float) : The area that produces the target MW
        """  

        func = lambda x: np.abs(self.mutate(A_helio=x).perform_calc(Ib=Ib, transmittance=transmittance, nG=nG).PowerMax_MW - target_MW)
        x0 = np.array(target_MW) # initial value 
        res = minimize(func, x0, method='nelder-mead',
                    options={'xatol': 1e-4, 'disp': False})
        return res.x[0]


    def financial_assessment(self, 
            oTow:OutputContainer,
            ee:EconomicEnvironment,
            csp_area_costs:float,
            power_block_cost:float,
            csp_energy_price:float,
            csp_discount_rate:float,
            lifetime:int|list=range(31)):
        """This function performs an economic analysis on the performance output of a csp

        Args:
            oTow (OutputContainer): _description_
            csp_area_costs (float): land costs in Euro/m2
            power_block_cost (float): power block cost in Euro/MW
            csp_energy_price (float): energy selling price in Euro/MWh
            csp_discount_rate (float): economic discount rate
            
            lifetime (int|list, optional): lifetime of the system. Defaults to range(30).
            
        Returns:
            _type_: _description_
        """
        if isinstance(lifetime, int):
            lifetime = list(range(lifetime+1))

        A_helio = self.A_helio_m2
        
        capital_csp_tow =A_helio* csp_area_costs + oTow.PowerMax_MW*power_block_cost
        revenue_csp_tow = cashflow(oTow.Energy_MWh, csp_energy_price, 
                                   fuel_energy=ee._Eoil, eff=0.4, fuel_price=-ee.oil_price,
                                   capital=capital_csp_tow)
                          
        cash_flow_tow = [-capital_csp_tow] + [revenue_csp_tow for i in lifetime[:-1]]
        
        # investment metrics
        dpb_tow = discounted_payback_period(csp_discount_rate, cash_flow_tow)
        npv_csp_tow = npf.npv (csp_discount_rate, cash_flow_tow)
        irr_csp_tow = npf.irr(cash_flow_tow )
        
        df = pd.DataFrame({'year': lifetime, 
                           'cash_flow': cash_flow_tow,
                           'discounted_cash_flow': [cf/(1+csp_discount_rate)**i for i,cf in enumerate(cash_flow_tow)]
                           }
                          )
        df['cumulative_cash_flow'] = df['discounted_cash_flow'].cumsum()
        
        return {
            'system_type': self._system_type,	
            'cash_flow_df': df,
            'scenario_params': oTow.scenario_params,
            'scenario_financial': {
                'PowerMax_MW': oTow.PowerMax_MW,
                'Energy_MWh': oTow.Energy_MWh,
                'CF': oTow.CF,
                'discounted_payback_period': dpb_tow, 
                'npv': npv_csp_tow,
                'irr': irr_csp_tow
                }
            }
        

# %% Old obsolete code
# def solarII(Ib:pd.Series,Trans:float,IAM:np.array,A_helio:float,Ar:float, 
#             nG:float = 0.97)->pd.Series:
#     """Calculates the power of the solar tower with heliostat
 
#     R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
#     Energy Conversion and Management. 60 (2012) 226–232.
#     J.E. Pacheco, R.W. Bradshaw, D.B. Dawson, W.D. la Rosa, R. Gilbert, S.H. Goods, P. Jacobs, M.J. Hale, S.A. Jones, G.J. Kolb, M.R. Prairie, H.E. Reilly, S.K. Showalter, L.L. VANT-HULL, 
#     Final Test and Evaluation Results from the Solar Two Project, n.d. https://core.ac.uk/reader/193342950 (accessed September 8, 2020).
 


#     Args:
#         Ib (pd.Series): direct irradiance
#         Trans (float): transmissivity 
#         IAM (np.array): incidence angle modifier
#         A_helio (float): heliostat area in m^2
#         Ar (float): receiver area in m^2
#         nG (float):  efficiency of generator Mosleh19 (Defaults to 0.97)

#     Returns:
#         pd.Series: power in MW
#     """
#     Effopt=100 # heliostat optical effiency [%] 65% pp.24 in Pacheco
#     #A_helio = 71140 + 10260 # total heliostat area [m2] pp.22 in Pacheco
#     R = 1 # reflectivity [%] 1 if IAM is IAM_tow(hoy)
#     Qin = Ib * R * Trans * IAM * A_helio * Effopt/100 # Eq. 17  in McGovern12
    
#     epsilon = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
#     sigma = 5.67 * 1e-8 # Stefan – Boltzman constant [W/m2K4]
#     Trec = 565 # the working fluid temperature in the receiver [oC] 565 oC 838K
#     Ta = 15 # ambient temperature close to the receiver [oC] 15 oC 288K
#     Tin = 290 # working fluid inlet temperature to the receiver [oC] 290 oC 563K
#     alpha = 1 # absorptivity of the receiver
#     #Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
    
#     Qrad = epsilon * sigma * Ar * (CtoK(Trec)**4-CtoK(Ta)**4)
#     hconv = CtoK(Trec)/60 + 5/3 # [W/m2K] convection coefficient Eq. 20 in McGovern12
#     '''D.L. Siebers, J.S. Kraabel, Estimating convective energy losses from solar central 
#     receivers, Sandia National Lab. (SNL-CA), Livermore, CA (United States), 1984. 
#     https://doi.org/10.2172/6906848.'''
#     Qconv = hconv * Ar * (CtoK(Trec) - CtoK(Ta))
    
#     Qnet = alpha * Qin - Qrad - Qconv # Eq. 8 in McGovern12
    
#     nR = 0.412 # SAM default

#     if Qnet.all() <= 0: #<<<<<<<<<<<<<<<<<<< check with Qin
#         P = 0
#     else:
#         P = Qnet * nR * nG
#     return P/1e6 # convert W to MW

#%% Alternative Incidence angle methods for towers


# def IAM_tow2(hoy:np.array=HOYS_DEFAULT) ->np.array : # polynomial fit, see file IAM.py for data
#     """Incidence angle modifier of Tower - elevation

#     Args:
#         hoy (np.array): hour of year

#     Returns:
#         np.array : Incidence angle modifier of Tower - elevation in rad
#     """    
#     return 1.66741484e-1 + 1.41517577e-2 * np.degrees(sgh.ele(hoy)) - 9.51787164e-5 * np.degrees(sgh.ele((hoy)))**2
