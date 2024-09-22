
import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from CombiCSP import HOYS_DEFAULT

class OutputContainer():
    """this is a container for the data and reshaping them. 

    **IMPORTANT NOTE**:    
    For simplicity reasons this assumes that the time index will always be hourly the days of the year

    #TODO consider giving (reference) to the original object that spawned this object. 
    # Alternatively use a dictionary for the parameters of the object

    """    
    hoy = HOYS_DEFAULT
    _df:pd.DataFrame = None
    scenario_params:dict = None
    system_type:str = None
    
    def __init__(self, power_df:pd.DataFrame, scenario_params:dict, system_type:str
                # ,	   hoys:np.ndarray = None
                 ):
        """constructor for the OutputContainer

        Args:
            data (np.array): time series
            A_helio (_type_): The (solar tower) heliostats area in m2
            Ctow (_type_): the ratio of area of heliostat to solar tower receiver 
        """
        # assert data.shape == hoys.shape, "Data and hoys should have the same shape"
        if not power_df.shape[0] == 8760:
            raise ValueError("Wrong time series")
        self.data_df = power_df
        self.scenario_params = scenario_params
        self.system_type = system_type
        assert isinstance(power_df, pd.DataFrame), "The data should be a pandas DataFrame"

        
    # def hour_power_arr(self):
    #     # raise PendingDeprecationWarning("This method will be removed in the future")
    #     return np.vstack((self.hoy, self.data_df))
    
    def power_data_as_array(self):
        """returns the power data as a numpy array

        Returns:
            np.array: the power data
        """        
        return self.data_df['Power_MW'].values
    
    def data4surf(self):
        """stacks the data in days and hours 

        Returns:
            _type_: _description_
        """        
        #TODO check if its possible to create a array of ???x24
        return np.vstack(self.power_data_as_array()).reshape((365,24))
    
    @property
    def PowerMax_MW(self)->float:
        """returns the maximum power of the Tower

        Returns:
            _type_: _description_
        """        
        return np.amax(self.power_data_as_array()) 

    @property
    def Energy_MWh(self)->float:
        """retunrs the total energy yield of the solar tower. 

        Returns:
           float : the total power per year in [MWh]
        """        
        return integrate.trapezoid(self.power_data_as_array()).round(2)
    
    @property
    def CF(self)->float: 
        """TODO Capacity factor?

        Returns:
            float: _description_
        """        
        return self.Energy_MWh / (8760 * self.PowerMax_MW)
    
    def as_df(self)->pd.DataFrame: 
        return self.data_df

    def summary_tower_data(self):    
        # TODO change the name (this does not apply to Tower data only)
        return self.scenario_params

    def print_summary(self):
        if self.system_type == 'tower':
            self._print_summary_tower()
        if self.system_type == 'trough':
            self._print_summary_trough()
        else:   
            print("No summary available")
    
    def _print_summary_tower(self):    
        
        str_out = ">>>>> Summary of results <<<<<<<<<\n" 
        str_out +="Reflective Solar Area : {} [m^2]\n".format(self.scenario_params.get('A_helio_m2', None))
        str_out +="ratio Reclective/Coll : {:.8g} \n".format( self.scenario_params.get('Ctow', None)) 
        str_out +="Max Power             : {:.2f} [MW]\n".format(self.PowerMax_MW)
        str_out +="Energy Yield          : {:.2f} [MWh]\n".format(self.Energy_MWh)
        str_out +="Concentration Factor  : {:.5g}  ".format(self.CF)
        print (str_out)

    def _print_summary_trough(self):    
        
        raise NotImplementedError("This method is not implemented yet")