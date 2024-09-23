# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
    @Credit: original functions based on G. Arnaoutakis
"""
#%%
import pathlib
import numpy as np
import pandas as pd
import pvlib

HOYS_DEFAULT = np.arange(1, 8761, 1) # hours of year

class SolarSystemLocation:
    def __init__(self, lat:float,  lon:float, 
                 dt_gmt_hr:float, 
                #  mer:float=0, # removed because it was not used. Probaly this refers to Local Standard Time Meridian (LSTM)
                 alt:float=0):
        """class that contains the location of the system 

        Args:
            lat (float): latitude of system in degrees (E-W in [-180,180])
            lon (float): longitude of system in degrees (N-S in [-90,90])
            mer (float):  for Greece check to replace with 15 * dt_gmt (TODO better description) see https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time
            dt_gmt_hr (float): time difference between Greenwich Mean Time in hours
            alt (float): altitude of system (in m) (default: 0  - sea level)
        """        
        self.dt_gmt = dt_gmt_hr # 
        self.lat_deg = lat # Crete
        self.lon_deg = lon # Crete 35.2401° N, 24.8093° E [east negative, west positive
        # self.mer = mer # TODO REMOVE mer from the init 
        # self.alt = alt # altitude of the system in meters above sea level
        

    @property
    def lat_rad(self)->float:
        """return latitude in radians

        Returns:
            _type_: _description_
        """        
        return np.radians(self.lat_deg)
    
    @property	
    def long_rad(self)->float:
        """return longitude in radians

        Returns:
            float: _description_
        """        
        return np.radians(self.lon_deg)
    
    def air_mass(self, hoy:np.array=HOYS_DEFAULT, method:str= 'wiki' ):
        """wrapper function for the different air mass
        
        Args:
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
            method (str, optional): _description_. Defaults to 'wiki'.

        Returns:
            _type_: _description_
        """    
        
        dic = {'wiki': self._AM_wiki,
            'Kasten': self._AM_Kasten,
            'Kasten-Young': self._AM3_KastenYoung,
            'Schoenberg' : self._AM_Shoenberg
            }
        if  method not in dic.keys():
            raise(ValueError(f"method can be [{dic.keys()}]"))
        return dic.get(method,None)(hoy)

    def _AM_wiki(self, hoy:np.array=HOYS_DEFAULT): # Air mass https://en.wikipedia.org/wiki/Air_mass_(solar_energy)
        AM =  1 / np.cos(self.z_rad(hoy))
        return AM
    
    def _AM_Kasten(self, hoy:np.array=HOYS_DEFAULT):
        '''F. Kasten, A new table and approximation formula for the relative optical air mass, 
        Arch. Met. Geoph. Biokl. B. 14 (1965) 206–223. https://doi.org/10.1007/BF02248840.'''
        return 1 / (np.cos(self.z_rad(hoy)) + 0.6556 * (6.379 - self.z_rad(hoy))**-1.757)
    
    def _AM3_KastenYoung(self, hoy:np.array=HOYS_DEFAULT):
        '''F. Kasten, A.T. Young, Revised optical air mass tables and approximation formula, 
        Appl. Opt., AO. 28 (1989) 4735–4738. https://doi.org/10.1364/AO.28.004735.'''
        return 1 / (np.cos(self.z_rad(hoy)) + 0.50572 * (6.07995 - self.z_rad(hoy))**-1.6364)
    
    def _AM_Shoenberg(self, hoy:np.array=HOYS_DEFAULT):
        '''E. Schoenberg, Theoretische Photometrie, in: K.F. Bottlinger, A. Brill, E. Schoenberg, 
        H. Rosenberg (Eds.), Grundlagen der Astrophysik, Springer, Berlin, Heidelberg, 1929: pp. 1–280. 
        https://doi.org/10.1007/978-3-642-90703-6_1.'''
        Re = 6371 # radius of the Earth [in km]
        yatm = 9 # effective height of the atmosphere [in km]
        r = Re / yatm
        return np.sqrt((r * np.cos(self.z_rad(hoy)))**2 + 2 * r + 1) - r * np.cos(self.z_rad(hoy))

    def tsol(self, hoy:np.array=HOYS_DEFAULT): # solar time [in decimal hours] introduce if function for east/west<<<<<<<<<
        """returns solar time 

        CALCULATION NEGLECTS DAYLIGHT SAVING ON SUMMER
        https://www.pveducation.org/pvcdrom/properties-of-sunlight/the-suns-position
        hoy + ((lat-mer)/15 + EoT(hoy))/60 # [60 min/h]
        pp.5, Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
        https://doi.org/10.1201/b22190

        Args:
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            _type_: _description_
        """
        time_correction_factor = 4*(self.lon_deg-15*self.dt_gmt) + EoT(hoy)
        return hoy + (time_correction_factor)/60 # [60 min/h]


    def W(self, hoy:np.array=HOYS_DEFAULT): # solar hour angle [in degrees]
        #TODO rename to Omega
        ''' Returns the solar hour in degrees. (the angle between the sun and the local meridian.)

        for tsol = 12h it should be ω = 0ο and for the solar time range
        tsol = 0 – 24h the solar hour angle ranges from 0o to ±180
        '''
        solar_time = self.tsol(hoy)
        return 15 * (solar_time - 12) # 360deg/24h = 15deg/h

    def ele_rad(self, hoy:np.array=HOYS_DEFAULT): # solar elevation angle or solar height [in radians]
        return np.arcsin(np.cos(np.radians(self.lat_deg)) * np.cos(delta_rad(hoy)) * np.cos(np.radians(self.W(hoy))) \
            + np.sin(np.radians(self.lat_deg)) * np.sin(delta_rad(hoy)))

    def z_rad(self, hoy:np.array=HOYS_DEFAULT):
        """Returns the solar zenith angle in radians

        solar zenith angle [in radians] https://en.wikipedia.org/wiki/Solar_zenith_angle   
        Args:
            hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

        Returns:
        np.array: solar zenith angle in radians"""

        return np.arccos(np.cos(np.radians(self.lat_deg)) * np.cos(delta_rad(hoy)) * np.cos(np.radians(self.W(hoy))) 
        + np.sin(np.radians(self.lat_deg)) * np.sin(delta_rad(hoy)))

    def azim_rad(self, hoy:np.array=HOYS_DEFAULT)->np.array: 
        """Returns the solar azimuth angle in radians

        Args:
            hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

        Returns:
            np.array: solar azimuth angle in radians
        """   
        return np.arcsin(np.cos(delta_rad(hoy)) * np.sin(np.radians(self.W(hoy))) / np.cos(self.ele_rad(hoy)))

    def Ib_from_csv(self,FNAME:pathlib.Path)->pd.Series:
        """loads data from a local csv folder

        Args:
            FNAME (pathlib.Path): path to csv file
        """        
        pvgis_data = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
        self.Ib = pvgis_data.loc[:,'Gb(n)']
        return self.Ib

    def Ib_from_pvgis(self, tz: str = 'Europe/Athens', 
        altitude:float = 400, 
        name = 'Default')->pd.Series:
        """Connects to PVGIS database and downloads the file.

        Args:
            tz (str, optional): Timezone. Defaults to 'Europe/Athens'.
            altitude (float, optional): Site altitude in [m]. Defaults to 400.
            name (str, optional): Name of the site. Defaults to 'Default'.

        Returns:
            pd.DataFrame: A DataFrame containing the following columns, with hourly data as index:
                - ghi: Global horizontal irradiance [W/m^2].
                - dni: Direct normal irradiance [W/m^2].
                - dhi: Diffuse horizontal irradiance [W/m^2].

            Note: The times of data are for the year 2020.
            
        """

        times = pd.date_range(start='2020-01-01', periods=8760, freq='1H', tz=tz) #end='2020-12-31', 
        #TODO allow for modification of the times of start and end
        solpos = pvlib.solarposition.get_solarposition(times, self.lat_deg, self.lon_deg)
        apparent_zenith = solpos['apparent_zenith']
        airmass = pvlib.atmosphere.get_relative_airmass(apparent_zenith)
        pressure = pvlib.atmosphere.alt2pres(altitude)
        airmass = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)
        linke_turbidity = pvlib.clearsky.lookup_linke_turbidity(times,  self.lat_deg, self.lon_deg)
        dni_extra = pvlib.irradiance.get_extra_radiation(times)
        pvlib_data = pvlib.clearsky.ineichen(apparent_zenith, airmass, linke_turbidity, altitude, dni_extra)
        self.Ib= pvlib_data.dni
        return pvlib_data


#TODO: consider creating a system/Unit parameters
#  def thetai(hoy:np.array=HOYS_DEFAULT, inclination=90, azimuths=0): # incidence angle [in radians]
#
#     g = deg(azim(hoy)) - azimuths # if surface looks due S then azimuths=0
#     return np.arccos(np.cos(ele(hoy)) * np.sin(np.radians(inclination)) * np.cos(np.radians(g))
#         + np.sin(ele(hoy)) * np.cos(np.radians(inclination)))

# ssCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)


def get_pvgis_tmy_data(sysloc:SolarSystemLocation)->pd.DataFrame:
    """function that collects data from online. 

    #TODO this should be a class containing tmy data which can either be retrieved from a file on disk or PVGIS database
    

    Args:
        sysloc (_type_): _description_

    Returns:
        pd.DataFrame: _description_
    """    
    latitude = sysloc.lat_deg
    longitude = sysloc.long
    OUTPUTFORMAT = 'json'

    dat = pvlib.iotools.get_pvgis_tmy(latitude, longitude, outputformat=OUTPUTFORMAT , usehorizon=True, userhorizon=None, 
        startyear=None, endyear=None, url='https://re.jrc.ec.europa.eu/api/', map_variables=True, timeout=30)

    df = dat[0]
    # dat[1] # monts of year for data 
    # dat[2] # metadata for the data set
    # dat[3] # variable explanation
    return df

#%% Equations of time

def EoT(hoy:np.array=HOYS_DEFAULT): # equation of time [in minutes]
    """
    Equation of time from Duffie & Beckman and attributed to Spencer
    (1971) and Iqbal (1983) [in minutes]
    """
    gamma = 360*(hoy-1)/365 #FIXME  hoy is in hours, while the equation should be in days
    return 2.2918*(0.0075+0.1868*np.cos(np.radians(gamma))-3.2077*np.sin(np.radians(gamma)) \
        -1.4615*np.cos(np.radians(2*gamma))-4.089*np.sin(np.radians(2*gamma)))

# TODO this do not work due to dayofyear --- uncomment when this is clear.
# def _calculate_simple_day_angle(dayofyear, offset=1):
#     """simple method for calculating the solar angle

#     Args:
#         dayofyear (_type_): _description_
#         offset (int, optional): _description_. Defaults to 1.

#     Returns:
#         _type_: solar angle  in radians 
#     """    
#     return (2. * np.pi / 365.) * (dayofyear - offset)

# def EoTS(hoy): # Equation of time from Duffie & Beckman and attributed to Spencer
#     #(1971) and Iqbal (1983) [in minutes]
#     day_angle = _calculate_simple_day_angle(dayofyear)
#     return (1440.0 / 2 / np.pi) * (0.0000075 +
#     0.001868 * np.cos(day_angle) - 0.032077 * np.sin(day_angle) -
#     0.014615 * np.cos(2.0 * day_angle) - 0.040849 * np.sin(2.0 * day_angle))

# def EoTPVCDROM(hoy): # equation of time [in minutes]
#     # PVCDROM: http://www.pveducation.org/pvcdrom/2-properties-sunlight/solar-time
#     # Soteris A. Kalogirou, "Solar Energy Engineering Processes and
#     # Systems, 2nd Edition" Elselvier/Academic Press (2009).
#     bday = _calculate_simple_day_angle(dayofyear) - (2.0 * np.pi / 365.0) * 80.0
#     return 9.87 * np.sin(2.0 * bday) - 7.53 * np.cos(bday) - 1.5 * np.sin(bday)


#%% ===================================== earth declination angles
def eda(hoy:np.array=HOYS_DEFAULT, method:str= 'wiki' ):
    """wrapper function for the different earth declination angles functions

    Useful reading [sunpos.py](https://levelup.gitconnected.com/python-sun-position-for-solar-energy-and-research-7a4ead801777)

    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
        method (str, optional): _description_. Defaults to 'wiki'.

    Returns:
        _type_: _description_
    """    
    
    dic = {'wiki':delta_rad,
           'Katsaprakakis': d2,
           '-81': d3,
           'pveducation' :d1
           }
    if  method not in dic.keys():
        raise(ValueError(f"method can be [{dic.keys()}]"))
    return dic.get(method,None)(hoy)

def d1(hoy:np.array=HOYS_DEFAULT): # earth declination angle [in degrees]
    """earth declination angle 

    The +10 comes from the fact that the winter solstice occurs before the start
    of the year. The equation also assumes that the suns orbit is a perfect circle and 
    the factor of 360/365 converts the day number to a position in the orbit.

    Args:
        hoy (np.array, optional): an array in hours of year . Defaults to HOYS_DEFAULT.

    Returns:
        _type_: earth inclination angle in degrees
    """
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle#footnote1_osno74c'''
    return -23.45 * np.cos( np.radians(360*(hoy+10*24)/365*24))

def d2(hoy:np.array=HOYS_DEFAULT):
    '''pp.6, Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190.'''
    return 23.45 * np.sin(np.radians(360*(hoy+284*24)/365*24))

def d3(hoy:np.array=HOYS_DEFAULT):
    return 23.45 * np.sin(np.radians(360*(hoy-81*24)/365*24))

def delta_rad(hoy:np.array=HOYS_DEFAULT): 
    '''Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    ecliptic in radians but without limiting the angle to be less than 2*Pi
    (i.e., the result may be greater than 2*Pi)
    
    [in radians] https://en.wikipedia.org/wiki/Sunrise_equation
    '''
    dOmega = 2.1429 - 0.0010394594 * hoy
    dMeanLongitude = 4.8950630 + 0.017202791698 * hoy
    dMeanAnomaly = 6.2400600 + 0.0172019699 * hoy
    dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) 
    + 0.00034894 * np.sin( 2 * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * hoy
    + 0.0000396 * np.cos(dOmega)
    '''Calculate celestial coordinates ( right ascension and declination ) in radians
    but without limiting the angle to be less than 2*Pi (i.e., the result may be
    greater than 2*Pi)'''
    dSin_EclipticLongitude = np.sin( dEclipticLongitude )
    dY = np.cos( dEclipticObliquity ) * dSin_EclipticLongitude
    dX = np.cos( dEclipticLongitude )
    dRightAscension = np.arctan2( dY,dX )
    if dRightAscension.any() < 0: dRightAscension = dRightAscension + 2*np.pi
    return np.arcsin( np.sin( dEclipticObliquity ) * dSin_EclipticLongitude )



def ineichen(latitude:float = 35, longitude: float = 24, 
    tz: str = 'Europe/Athens', 
    altitude:float = 400, 
    name = 'Ierapetra' )->pd.DataFrame:
    """function that calls the pvgis and obtains the meteorological data

    Args:
        latitude (float, optional): _description_. Defaults to 35.
        longitude (float, optional): _description_. Defaults to 24.
        tz (str, optional): _description_. Defaults to 'Europe/Athens'.
        altitude (float, optional): _description_. Defaults to 400.
        name (str, optional): _description_. Defaults to 'Ierapetra'.

    Returns:
        pd.DataFrame: a Dataframe containing the following columns with index the hourly data
            ghi:  global horizontal irradiance[W/m^2]
            dni:  direct normal irradiance [W/m^2]
            dhi:  diffuse horizontal irradiance [W/m^2]
    """   
    times = pd.date_range(start='2020-01-01', periods=8760, freq='1H', tz=tz) #end='2020-12-31', 
    solpos = pvlib.solarposition.get_solarposition(times, latitude, longitude)
    apparent_zenith = solpos['apparent_zenith']
    airmass = pvlib.atmosphere.get_relative_airmass(apparent_zenith)
    pressure = pvlib.atmosphere.alt2pres(altitude)
    airmass = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)
    linke_turbidity = pvlib.clearsky.lookup_linke_turbidity(times, latitude, longitude)
    dni_extra = pvlib.irradiance.get_extra_radiation(times)
    return pvlib.clearsky.ineichen(apparent_zenith, airmass, linke_turbidity, altitude, dni_extra)