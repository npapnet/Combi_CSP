#%% [markdown]
# # Scope
# this module contains a class that is able to 
# - handle a system which contains a combination of many Solar Troughs and Tower
# - output a container similar in functionality with the Solar Trough and Tower. 
#%%
import logging 
from CombiCSP import SolarTowerCalcs, SolarTroughCalcs
#%% 
class CSPSystemCombination():
    """class that handles the calculation of a CSP system.
    """    
    def __init__(self, tow_lst:list= None, trough_lst:list= None):
        """constructor 

        Args:
            tow_lst (list): list of solar tower objects
            trough_lst (list): list of solar trough objects
        """
        self._tow_lst = tow_lst if tow_lst is not None else tow_lst
        self._trough_lst = trough_lst if trough_lst  is not None else trough_lst 
    
    def add_tower(self, tow_obj):
        """_summary_

        Args:
            tow_obj (_type_): _description_
        """
        assert isinstance(tow_obj, SolarTowerCalcs)
        self._tow_lst.append(tow_obj)

    def add_trough(self, trough_obj):
        """_summary_

        Args:
            trough_obj (_type_): TroughCalcs obj
        """
        assert isinstance(trough_obj, SolarTroughCalcs)
        self._trough_lst.append(trough_obj)

    def remove_trough(self, position_index:int):
        """function to remove trough

        Args:
            position_index (_type_): positional index
        """
        try:
            self._trough_lst.pop(position_index)
        except:
            logging.error('Could not remove object from trough list')        
        

    def remove_tower(self, position_index:int):
        """function to remove trough

        Args:
            position_index (_type_): positional index
        """
        try:
            self._tow_lst.pop(position_index)
        except:
            logging.error('Could not remove object from tower list')        
        

    def perform_calc(self, hoy, Ib, tow_args:dict, trough_args:dict):
        """_summary_

        Args:
            hoy (_type_): _description_
            Ib (_type_): _description_
            tow_args (dict): _description_
            trough_args (dict): _description_

        Raises:
            Exception: _description_
        """        
        raise Exception('Not implemented yet')
        res =0
        return res
    