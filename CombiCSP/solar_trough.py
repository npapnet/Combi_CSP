# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis, G. Arnaoutakis
    @Date: 2022/07/02
    @Credit: original functions from G. Arnaoutakis
"""
#%%
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import numpy_financial as npf

from CombiCSP import OutputContainer, CtoK, HOYS_DEFAULT
from CombiCSP.solar_system_location import SolarSystemLocation, delta_rad

from CombiCSP.financial.economics import cashflow, discounted_payback_period, EconomicEnvironment


STEFAN_BOLTZZMAN_CONSTANT = 5.67 * 1e-8 # [W/m2K4] Stefan – Boltzman constant
class SolarTroughCalcs():
    _system_type:str = 'trough'
    _system_alignment:str = None # NS or EW (currently not implemented)

    _hourly_results : OutputContainer = None
    N:int = None # number of units
    L:float= None # solar trough unit length
    Ws:float = None # width between rows
    receiver_dia_m:float = None # tube outer diameter
    collector_width_m:float = None # collector width
    foc_len:float = None # focal length

    # calculation constants
    optical_efficiency:float = 75 # [%] Optical efficiency of collector 74% INDITEP in pp.7 Fraidenraich13, pp.4 Zarza06
    epsilon:float = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
    #Tr = 350 # [oC] the working fluid temperature in the receiver DISS pp.3,7 in Zarza04
    alpha:float = 1 # absorptivity of the receiver
    nR:float = 0.375 # isentropic efficiency Salazar17

    alignment_calculator = None # Class that calculates the power of the system based on alignment EW or NS

    def __init__(self,
        foc_len = 0.88 
        ,N = 1800 
        ,L = 25 
        ,Ws = 18 
        ,Wr = 0.07 
        ,Wc = 5.76
        , slobj:SolarSystemLocation =  None
        ):
        """_summary_

        Args:
            foc_len (float): [m] focal length CSPP T.1 in Mosleh19. Defaults to 0.88.
            N (int):  # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine. Defaults to 1800.
            L (int): [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  . Defaults to 25.
            Ws (int): [m] width between rows 18 INDITEP in pp.6 Fraidenraich13,  pp.5 Zarza06. Defaults to 18.
            Wr (float): tube outer diameter [m]. Defaults to 0.07.
            Wc (float): collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM. Defaults to 5.76.
        """        
        if isinstance(slobj, SolarSystemLocation):
            self._sl = slobj
        else:
            raise ValueError('Site Location has not been provided')

        self.foc_len = foc_len
        self.N = N
        self.L = L
        self.Ws = Ws
        self.receiver_dia_m = Wr
        self.collector_width_m = Wc 

    def print_system_summary(self):
        print(f"""
--------  System Parameters
-> Focal Lentth     (foc_len) = {self.foc_len}  
-> Number of units  (   N   ) = {self.N}      
-> Unit length      (   L   ) =  {self.L} 
-> Unit spacing?    (   Ws  ) = {self.Ws} 
-> Receiver width   (   Wr  ) = {self.receiver_dia_m}
-> Collector width  (   Wc  ) = {self.collector_width_m}
--------  Derived Quantities
-> Collector area   (   Ac  ) = {self.Ac()}
-> Receiver area    (   Ar  ) = {self.Ar()}
-> Concetration f   (   Cg  ) = {self.Cg}
        """)        

    @property
    def area(self):
        """returns the collector area. 

        update Wc for geometry, see A. Rabl, Comparison of solar concentrators, Solar Energy. 18 (1976) 93–111.
        
        Args:
            Wc (float): width of solar collector in m^2
            L (float): length of solar collector in m^2
            N (int): quantity of solar collectors. 

        Returns:
            _type_: total collector area.


        Returns:
            float: The collector area in [m^2]
        """        
        return self.collector_width_m * self.L * self.N
    
    def Ac(self)->float:
        """collector area in m^2

        (Assumption): it is the width*length times the number of units

        Returns:
            float: collector area in m^2
        """        
        return self.collector_width_m * self.L * self.N

    def Ar(self)->float:
        """returns the receiver area. 

        Args:
            Wr (float): width of solar receiver in [m]
            L (float): length of solar receiver in [m]
            N (int): quantity of solar receivers. 

        Returns:
            float: total receiver area in m^2
        """    
        return self.receiver_dia_m * self.L * self.N
    
    @property
    def Cg(self)->float:
        """Geometrical Concentration (collector area/receiver area)

        Affected by the following properties:
            Wc (float): width of collector in  [m]
            Wr (float): width receiver in  [m]
            L (float): length in  [m]
            N (int): quantity of units []

        Returns:
            float: geometrical concentration of parabolic troughs [dimensionless]
        """

        return self.Ac() / self.Ar()

    def set_alignment(self, alignment:str):
        assert alignment in ['EW', 'NS'], 'Alignement should be one of ["EW", "NS"]'
        self._system_alignment = alignment
        if alignment == 'NS':
            self.alignment_calculator = SolarTroughNS(self)
        elif alignment == 'EW':
            self.alignment_calculator = SolarTroughEW(self) 
        else:
            raise Exception('Alignement should be one of ["EW", "NS"]')
        
    def perform_calc(self, alignment:str, Ib, Tr=318, hoy=HOYS_DEFAULT)->OutputContainer:
        self.set_alignment(alignment)
        return self.alignment_calculator.perform_calcs(Ib=Ib, Tr=Tr, hoy=hoy)

    def params_as_dict(self):
        dic = {
            'foc_len_m':self.foc_len,
            'no_units':self.N,
            'unit_length_m':self.L,
            'Ws':self.Ws,
            'receiver_dia_m':self.receiver_dia_m,
            'collector_width_m':self.collector_width_m,
            'area_m2':self.area,
            'Ac_m2':self.Ac(),
            'Ar_m2':self.Ar(),
            'Cg':self.Cg,
            'optical_efficiency':self.optical_efficiency,
            'epsilon':self.epsilon,
            'alpha':self.alpha,
            'nR':self.nR
        }
        return dic



    def thetai(self, hoy:np.array=HOYS_DEFAULT, inclination=90., azimuths=0.)->np.array: #
        """ Calculates the incidence angle [in radians]

        #TODO check whether there is a dependence between azimuth and NS and EW type of CSP

        Args:
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
            inclination (float, optional): _description_. Defaults to 90.
            azimuths (float, optional): _description_. Defaults to 0.

        Returns:
            np.array: _description_
        """        
        g = np.degrees(self._sl.azim_rad(hoy)) - azimuths # if surface looks due S then azimuths=0
        elev = self._sl.ele_rad(hoy)
        return np.arccos(np.cos(elev) * np.sin(np.radians(inclination)) * np.cos(np.radians(g)) 
            + np.sin(elev) * np.cos(np.radians(inclination)))


    def IAM_tro(self, hoy:np.array=HOYS_DEFAULT)->np.array: 
        """Incidence angle modifier of parabolic trough - equation1
        
        G.A. Salazar, N. Fraidenraich, C.A.A. de Oliveira, O. de Castro Vilela, M. Hongn, J.M. Gordon, 
        Analytic modeling of parabolic trough solar thermal power plants, Energy. 138 (2017) 1148–1156. 
        https://doi.org/10.1016/j.energy.2017.07.110.

        #TODO there are 4 different function for IAM_tro. They need to be consolidated in a single one and selected as an option.

        # thetai in radians

        Args:
            hoy (np.array): hour of year

        Returns:
            _type_: _description_
        """
        #TODO needs rad despite thetai(hoy) already in rad???
        thetas = self.thetai(hoy)
        return np.cos(np.radians(thetas)) + 0.02012 * thetas - 0.01030 * thetas**2

    
    def incident_energy_on_system(self,  Ib:pd.Series, hoy:np.array = HOYS_DEFAULT, alignment:str=None)->pd.Series:
        """Function that returns the incident energy on the system depending on the alignment

        Args:
            Ib (pd.Series): _description_
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
            alignment (str): can be 'NS' or 'EW'. Defaults to None (no changes in aligment).	

        Raises:
            Exception: _description_

        Returns:
            pd.Series: incident angle at different hour of the year. 
        """
        if alignment is not None:
            self.set_alignment(alignment)
        
        cos_thetas = self.alignment_calculator.costhetai(hoy)
        return Ib*self.IAM_tro(hoy) *cos_thetas

    def di_sst(self, hoy, Ib, costhetai, Tr, nG:float = 0.97
               , Ta = 15 # [oC] ambient temperature close to the receiver 15 oC 288K
            ,Tin = 200 # [oC] working fluid inlet temperature to the receiver DISS pp.3,7 in Zarza04
            )->pd.Series:
        """Calculates the total power of the parabolic system

        R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
        Energy Conversion and Management. 60 (2012) 226–232.
        E. Zarza, L. Valenzuela, J. León, K. Hennecke, M. Eck, H.-D. Weyers, M. Eickhoff, 
        Direct steam generation in parabolic troughs: Final results and conclusions of the DISS project, 
        Energy. 29 (2004) 635–644. https://doi.org/10.1016/S0360-5442(03)00172-5.


        Args:
            hoy (np.array): hour of year
            Ib (np.array): hour of year
            costhetai (_type_): cosine function [rad] (not necessary)
            IAM (_type_): incidence angle modifier
            Tr (float): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04 
            nG (float): Generator efficiency [dimensionless]

        Returns:
            _type_: power in [MW].
        """    
        
        incident_energy = self.incident_energy_on_system(Ib=Ib, hoy=hoy) #( this is based on alignment)
        # Thermal Energy input (in the duration of an hour)
        # Qin = Ib * costhetai* IAM * self.Ac() * self.optical_efficiency/100 # Eq. 4  in McGovern12
        Qin = incident_energy * self.Ac() * self.optical_efficiency/100 # Eq. 4  in McGovern12
        
        # Radiation losses
        Qrad = self.epsilon * STEFAN_BOLTZZMAN_CONSTANT * self.Ar() * (CtoK(Tr)**4-CtoK(Ta)**4) # check model from Broesamle
        # Convection losses
        Qconv = 0
        # Energy Balance
        Qnet = self.alpha * Qin - Qrad - Qconv # Eq. 8 in McGovern12
        
        # ==   Carnot calulation for the power output of the system (not used)
        # Turbine 
        Tcond = 30 # [oC] condenser temperature
        Ts = 565 # [oC] steam temperature
        n_Carnot = 1 - (CtoK(Tcond)/CtoK(Ts)) # Eq. 15  in McGovern12
        # =============================

        # power output 
        if Qnet.all() <= 0:
            P = 0
        else:
            P = Qnet * self.nR * nG
        return P/1e6 # convert W to MW

    def calculate_nominal_power_MW(self, Tr:float=318, T_amb:float=15, nG:float= 0.97)->pd.Series:
        """
        Calculates the nominal power of the solar trough system

        The nominal power in calculated based on:
        - the normal irradiance Ib set to 1000 W/m2 
        - the working fluid temperature in the receiver Tr (Tr)
        - the ambient temperature T_amb (T_amb) set to 318
        - the generator efficiency nG (self.nG) set to 15 C
        - the optical efficiency of the collector (self.optical_efficiency)
        """
        Ib = 1000 # [W/m2] beam irradiance
        # Tr = 318 # [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04
        # T_amb = 15 # [oC] ambient temperature close to the receiver 15 oC 288K
        # hoy = HOYS_DEFAULT

        Q_in = self.Ac()*self.optical_efficiency/100 * Ib
        Qrad = self.epsilon * STEFAN_BOLTZZMAN_CONSTANT * self.Ar() * (CtoK(Tr)**4-CtoK(T_amb)**4) # check model from Broesamle    
        Q_conv = 0
        Q_net = self.alpha * Q_in - Qrad - Q_conv
        P_MW = Q_net * self.nR * nG /1e6
        return P_MW

 
        
    def mutate(self,foc_len = None 
        ,N = None 
        ,L = None 
        ,Ws = None
        ,Wr = None 
        ,Wc = None
        , slobj:SolarSystemLocation =  None)->'SolarTroughCalcs':
        """ produces a mutated solar trough object, with different properties (When set). 

        Args:
            foc_len (_type_, optional): _description_. Defaults to None.
            N (_type_, optional): _description_. Defaults to None.
            L (_type_, optional): _description_. Defaults to None.
            Ws (_type_, optional): _description_. Defaults to None.
            Wr (_type_, optional): _description_. Defaults to None.
            Wc (_type_, optional): _description_. Defaults to None.
            slobj (SolarSystemLocation, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """        
        foc_len = self.foc_len if foc_len is None else foc_len
        N = self.N if N is None else N
        L = self.L if L is None else L
        Ws = self.Ws if Ws is None else Ws
        Wr = self.receiver_dia_m if Wr is None else Wr
        Wc = self.collector_width_m if Wc is None else Wc
        slobj = self._sl if slobj is None else slobj
        return SolarTroughCalcs(
            foc_len = foc_len 
            ,N = N 
            ,L = L 
            ,Ws = Ws 
            ,Wr = Wr 
            ,Wc =Wc 
            , slobj=  slobj
                                )

    def find_no_units_for_nominal_power_MW(self,
        target_MW :float,
        Tr=318.)->int:
        """finds the required area based on the other parameters of the solar trough

        Ib: 1000 W/m2
        Tr: 318oC
        T_amb: 15oC
        nG: 0.97
        """
        func = lambda x: np.abs(self.mutate(N=x).calculate_nominal_power_MW(Tr=318, T_amb=15, nG= 0.97) - target_MW)
        x0 = np.array(self.N) # initial value
        res = minimize(func, x0, method='nelder-mead',
                    options={'xatol': 1e-4, 'disp': False})
        return np.ceil(res.x[0])
    
    
    def find_units_for_max_MW(self, 
        target_MW :float,
        alignment:str, 
        Ib:pd.Series, 
        Tr=318., hoy=HOYS_DEFAULT)->int:
        """finds the required number of units based on the other parameters of the solar trough

        Args:
            target_MW(float): max MW targed
            alignment(str): the alignment of the solar trough 
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (_type_, optional): _description_. Defaults to HOYS_DEFAULT.
        
        returns:
            (int) : The area that produces the target MW
        """

        func = lambda x: np.abs(self.mutate(N=x).perform_calc(alignment=alignment, Ib=Ib, Tr=Tr).PowerMax_MW - target_MW)
        x0 = np.array(self.N) # initial value 
        res = minimize(func, x0, method='nelder-mead',
                    options={'xatol': 1e-4, 'disp': False})
        return np.ceil(res.x[0])


    def financial_assessment(self, 
        oTr:OutputContainer,
        ee:EconomicEnvironment,
        csp_area_costs:float,
        power_block_cost:float,
        csp_energy_price:float,
        csp_discount_rate:float,
        lifetime=range(31)):
        """This function performs an economic analysis on the performance 
        output of a csp

        Args:
            oTr (OutputContainer): references the trough object
            ee (Economic_environment): economic environment 
            csp_area_costs (_type_): csp area costs
            power_block_cost (_type_): cost of the power block
            csp_energy_price (_type_): energy prices
            csp_discount_rate (_type_): discount rate
            lifetime (_type_, optional): _description_. Defaults to (0,1,2.. 30).

        Returns:
            _type_: _description_
        """
        if isinstance(lifetime, int):
            lifetime = list(range(lifetime+1))

        area = oTr.scenario_params.get('area_m2')

        capital_csp_tro = area * csp_area_costs + oTr.PowerMax_MW*power_block_cost
        revenue_csp_tro = cashflow(Ecsp=oTr.Energy_MWh,csp_price=csp_energy_price,
                                   fuel_energy=ee._Eoil,eff=0.4,fuel_price=-ee.oil_price,
                                   capital=capital_csp_tro)
        
        cash_flow_tro = [-capital_csp_tro] + [revenue_csp_tro for _ in lifetime[:-1]]

        # investment metrics
        dpb_tro = discounted_payback_period(csp_discount_rate, cash_flow_tro)
        npv_csp_tro = npf.npv(csp_discount_rate, cash_flow_tro)
        irr_csp_tro = npf.irr(cash_flow_tro)
        
        df = pd.DataFrame({'year': lifetime, 
                           'cash_flow': cash_flow_tro,
                           'discounted_cash_flow': [cf/(1+csp_discount_rate)**i for i,cf in enumerate(cash_flow_tro)]
                           }
                          )
        df['cumulative_cash_flow'] = df['discounted_cash_flow'].cumsum()

        return {
            'system_type': oTr.system_type,
            'cash_flow_df': df,
            'scenario_params': oTr.scenario_params,
            'scenario_financial': {
                'PowerMax_MW': oTr.PowerMax_MW,
                'Energy_MWh': oTr.Energy_MWh,
                'CF': oTr.CF,
                'discounted_payback_period': dpb_tro, 
                'npv': npv_csp_tro,
                'irr': irr_csp_tro}
        }



#region NS Trough  
class SolarTroughNS():
    _system_type:str = 'trough_NS'

    def __init__(self, Trough:SolarTroughCalcs):
        self.Trough = Trough
        self._sl = self.Trough._sl
    
    def perform_calcs(self, Ib, Tr=318., hoy=HOYS_DEFAULT)->OutputContainer:
        """Calculation for a solar trough oriented NS for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """        
        system_type = self._system_type
        lat_rad = self._sl.lat_rad
        #Parabolic trough cosine function in North-South orientation
        #   Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
        #   Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.
 
        power_data = self.Trough.di_sst(hoy=hoy, Ib=Ib,costhetai=self.costhetai(hoy),
                      Tr=Tr)
        df = pd.DataFrame({'HOY':hoy,'Ib_n':Ib, 'Power_MW':power_data})
        scenario_params = {'system_type': system_type}
        scenario_params.update( self.Trough.params_as_dict())

        output  = OutputContainer(power_df = df, scenario_params=scenario_params, system_type=system_type)
        return output

    def costhetai(self, hoy)->np.array:
        lat_rad = self._sl.lat_rad
        deltas = delta_rad(hoy)
        ws_rad = np.radians(self._sl.W(hoy))
        return np.cos(deltas) * (np.sin(ws_rad)**2 + 
                (np.cos(lat_rad) *  np.cos(ws_rad) + np.tan(deltas) * np.sin(lat_rad))**2)**0.5
#endregion

#region EW Trough  
class SolarTroughEW():
    _system_type:str = 'trough_EW'

    def __init__(self, Trough:SolarTroughCalcs):
        self.Trough = Trough
        self._sl = self.Trough._sl
    

    def perform_calcs(self, Ib, Tr=318, hoy=HOYS_DEFAULT)->OutputContainer:
        """Calculation for a solar trough oriented EW for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """ 
        system_type= self._system_type 
        # IAM = self.Trough.IAM_tro(hoy)
        
        #Parabolic trough cosine function in East West orientation
        #    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
        #   Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.
        
        # costhetai_EW_arr = self.costhetai(hoy)
        power_data = self.Trough.di_sst(hoy = hoy, Ib=Ib,
                                        costhetai= self.costhetai(hoy), 
                                        Tr=Tr)
        df = pd.DataFrame({'HOY':hoy,'Ib_n':Ib, 'Power_MW':power_data})
        scenario_params = {'system_type': system_type}
        scenario_params.update( self.Trough.params_as_dict())
        output  = OutputContainer(power_df = df, scenario_params=scenario_params, system_type=system_type)
        return output

    def costhetai(self, hoy):
        deltas = delta_rad(hoy)
        ws_rad = np.radians(self._sl.W(hoy))
        return  np.cos( deltas) * (np.cos(ws_rad)**2 + np.tan(deltas**2))**0.5
    

#endregion



#%% 
#region obsolete code
# Incidence angle methods for troughs
# def IAM_tro2(hoy:np.array=HOYS_DEFAULT):
#     '''N. Fraidenraich, C. Oliveira, A.F. Vieira da Cunha, J.M. Gordon, O.C. Vilela, 
#     Analytical modeling of direct steam generation solar power plants, Solar Energy. 98 (2013) 511–522. 
#     https://doi.org/10.1016/j.solener.2013.09.037.
#     citing M.G. B, E.L. F, R.O. A, A.E. A, W.S. C, A.S. C, E.Z. E, P.N. B, 
#     EUROTROUGH- Parabolic Trough Collector Developed for Cost Efficient Solar Power Generation, in: n.d.'''
#     #return 1 - 0.00044 * thetai(hoy) / cos(thetai(hoy)) - 0.00003 * (thetai(hoy))**2 / cos(thetai(hoy)) # needs rad despite thetai(hoy) already in rad???
#     '''De Luca15
#     citing M.J. Montes, A. Abánades, J.M. Martínez-Val, M. Valdés, 
#     Solar multiple optimization for a solar-only thermal power plant, using oil as heat transfer fluid in the parabolic trough collectors, 
#     Solar Energy. 83 (2009) 2165–2176. https://doi.org/10.1016/j.solener.2009.08.010.
#     citing B, M.G.; F, E.L.; A, R.O.; A, A.E.; C, W.S.; C, A.S.; E, E.Z.; B, P.N. 
#     EUROTROUGH- Parabolic Trough Collector Developed for Cost Efficient Solar Power Generation.
#     '''
#     return np.cos(np.radians(thetai(hoy))) - 5.25097e-4 * thetai(hoy) - 2.859621e-5 * thetai(hoy)**2

# def IAM_tro3(hoy:np.array=HOYS_DEFAULT):
#     '''(Dudley, 1994)'''
#     return np.cos(np.radians(thetai(hoy))) - 0.0003512 * thetai(hoy) - 0.00003137 * (thetai(hoy))**2# thetai in degrees
#     '''pp.26 A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
#     University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
#     #return cos(rad(thetai(hoy))) + 0.000884 * thetai(hoy) - 0.00005369 * (thetai(hoy))**2 # needs rad despite thetai(hoy) already in rad???

# def IAM_tro4(hoy,foc_len,area,L):
#     '''eq.1 in H.J. Mosleh, R. Ahmadi, Linear parabolic trough solar power plant assisted with latent thermal energy storage system: 
#     A dynamic simulation, Applied Thermal Engineering. 161 (2019) 114204. https://doi.org/10.1016/j.applthermaleng.2019.114204.
#     J.A. Duffie, W.A. Beckman, Solar Engineering of Thermal Processes, Wiley, 1991.'''
#     return 1 - foc_len / L *(1 + area**2 / 48 * foc_len**2) * np.tan(np.rad(thetai(hoy))) # needs rad despite thetai(hoy) already in rad???

#%% Functions related to heat load of Conscentrated Soalr Collectors 
# # ============================================================================================
# def CSCUL(hoys): # CSC heat loss ex.4.2
#     """Concentrated Solar Collector heat loss ex.4.2

#     Returns:
#         _type_: _description_
#     """    
#     Dr = 0.07 # the receiver’s outer diameter [m]
#     Tr = 573 # the working fluid temperature in the receiver [K]
#     epsilonr = 0.25 # the emissivity of the receiver’s material
#     Dco = 0.1 # the outer diameter of the receiver’s cover [m]
#     t = 5*1e-3 # the receiver’s thickness [m]
#     u = 4 # wind velocity [m/s]
#     Ta = 288 # ambient temperature close to the receiver [K]
#     Tsky = 278 # the ambient temperature far from the receiver [K]
#     Ucond = 0.022 # the thermal transmittance factor for the heat transfer through conductivity [W/mK]
#     epsilonc = 0.90 # the emissivity of the cover’s material
#     kc = 1.45 # the receiver’s thermal conductivity factor [W/mK]
#     v = 1.456 * 1e-5 # kinematic viscosity of air with ambient temperature 15 οC [m2/s]
#     Tco = 306.36 # the temperature of the cover’s outer surface [K]
#     sigma = 5.67 * 1e-8 # Stefan – Boltzman constant [W/m2K4]
#     keff = 0 # the vacuum’s thermal conductivity factor [W/mK]
#     L = 25 # the length of the receiver [m]
#     Dci = Dco - (2 * t)
#     Re = u * Dco / v
    
#     Nu = 0.30 * Re**0.60
    
#     hw = Nu * Ucond / Dco
    
#     Qloss = np.pi * Dco * hw * (Tco - Ta) 
#     + epsilonc * np.pi * Dco * sigma * (Tco**4 - Tsky**4) * L
#     Tci = Qloss / (2 * np.pi * kc * L) * np.log(Dco/Dci) + Tco
#     Qloss2 = (2 * np.pi * keff * (CtoK(Tr) - Tci) / np.log(Dci/Dr) + 
#     (np.pi * Dr * sigma * (CtoK(Tr)**4 - Tci**4)) / (1/epsilonr+((1-epsilonc)*Dr/Dci)/epsilonc)) * L
#     for x in hoys:
#         UL = Qloss / (np.pi * Dr * L * (CtoK(Tr) - Ta))
#     return UL # [W/m2K]

# def CSCP(Tfi, hoy:np.array= HOYS_DEFAULT, fname:str="example_data/tmy_35.015_25.755_2005_2020.csv"): 
#     # CSC thermal power ex.4.3
#     """Concentrated Solar collect

#     Args:
#         Tfi (_type_): _description_
#         hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
#         fname (str, optional): _description_. Defaults to "example_data/tmy_35.015_25.755_2005_2020.csv".

#     Returns:
#         _type_: _description_
#     """    
#     '''see also Dikmen, E., Ayaz, M., Ezen, H.H., Küçüksille, E.U., Şahin, A.Ş., 2014. 
#     Estimation and optimization of thermal performance of evacuated tube solar collector system. 
#     Heat Mass Transfer 50, 711–719. https://doi.org/10.1007/s00231-013-1282-0
#     '''
    
#     try:
#         pvgis_data = pd.read_csv(fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')
#     except:
#         #TODO this is an exception until T = CSCP(Tr) is removed from this file (use tests)
#         pvgis_data = pd.read_csv("examples/"+fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')

#     Ib = pvgis_data.loc[:,'Gb(n)']
#     #S = 550 # incident solar radiation [W/m2]
#     #UL = 4.50 # heat losses total thermal transmittance factor [W/m2Κ]
#     Wc = 5.76 # width of the concentrating collector [m]
#     L = 12 # length of the concentrating collector [m]
#     cp = 2500 # fluid's specific heat capacity [J/kgK]
#     m = 8 # flow rate [kg/s] Mosleh19, pp.9
#     #Tfi = Tr # fluid's inlet temperature [C]
#     hfi = 400 # thermal convection factor for the heat transfer from the receiver to the working fluid [W/m2K]
#     k = 16 # stainless steel thermal conductivity factor [W/mK]
#     t = 5 * 1e-3 # receiver's tube thickness [m]
#     Ta = 15 # ambient temperature close to the receiver [C]
#     Dro = 0.07 # receiver’s outer diameter [m]
#     Dco = 0.1 # outer diameter of the receiver’s cover [m]
#     Ti = 270 #  average temperature of the receiver [C]
    
#     Ar = np.pi * Dro * L # receiver’s inner area [m2]
#     Aa = (Wc - Dco) * L # collector’s effective area [m2]
#     Dri = Dro - (2 * t) #  receiver’s inner diameter [m]
#     F = 1/CSCUL(hoy) / (1/CSCUL(hoy) + Dro/(hfi*Dri) + Dro * np.log(Dro/Dri)/(2*k))
#     FR = m * cp * (1 - np.exp(-Ar * CSCUL(hoy) * F / (m * cp))) / (Ar * CSCUL(hoy)) # heat removal factor - review precision <<<
#     Qu = Aa * FR * (Ib - Ar * CSCUL(hoy) * (Ti - Ta)/Aa) # final thermal power production [W]
#     Tfo = Tfi + (Qu / (m * cp)) # fluid’s outlet temperature [C]
#     Tro =  Tfi + (Qu * ((1/(np.pi * Dri * L * hfi)) + np.log(Dro/Dri)/(2*np.pi*k*L))) # receiver’s outer surface temperature [C]
#     DT = Tro - Tfi
#     return Tfo#Qu/1000 # convert W to kW

#endregion