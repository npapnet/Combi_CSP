# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""
#%% miscellaneous functions
import pathlib
import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from CombiCSP import HOYS_DEFAULT

def CtoK(c: float):
    """convert Celsius to Kelvin

    Args:
        c (float|int): the temperature in degrees Celsius

    Returns:
        float: the temperature in Kelvin
    """    
    k = c + 273
    return k
  


def heatmap2d_sns(data, title:str= '', figsize:tuple=(15,8)):
    ''' This function allows larger displays for the heatmap compared to the imshow.

    As a drawback it depends of seaborn.
    '''
    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax = sns.heatmap(data,ax =ax,
                cbar_kws={'label': 'Power [MW]'})
    ax.set_title ('')
    ax.set_xlabel ('Day of year')
    ax.set_ylabel ('Hour')
    ax.set_title(title)

def heatmap2d(arr: np.ndarray):
    plt.imshow(arr, cmap='hot', interpolation='gaussian')
    plt.xlabel('Day')
    plt.ylabel('Hour')
    #plt.savefig('maps.png') # <<<<<<<<<<<check operation
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)
    plt.ylabel('MW')#W/m${^2}$
    plt.show()
# %%

def compare_date(date_str:str, tower, trough, Ib:pd.Series, save_to_file:bool = True):
    """compares the incident energy on the collector 

    Args:
        tower (SolarTowerCalcs): _description_
        trough (SolarTroughCalcs): _description_
        date_str (str, optional): _description_. Defaults to '2020-08-26'.
        save_to_file (bool, optional): _description_. Defaults to True.
    """    
    date = pd.to_datetime(date_str)
    plt.figure()
    plt.plot(Ib)
    plt.plot(trough.incident_energy_on_system(alignment='EW', Ib=Ib) , label='Solar Trough EW')
    plt.plot(trough.incident_energy_on_system(alignment='NS', Ib=Ib),  label = 'Solar Trough NS')
    plt.plot(tower.incident_energy_on_system(Ib), label = 'Solar Tower')

    plt.xlim(date,date + pd.DateOffset(days=1))
    plt.ylim(0,800)
    plt.xlabel('Date-Time')
    plt.ylabel('Irradiance ($W/m^2$)')
    plt.legend(('DNI','Trough NS','Trough EW','Tower','tower','trough'))
    plt.xticks(rotation=30)


    Ens = integrate.trapezoid(trough.incident_energy_on_system(alignment='NS', Ib=Ib)).round(2)
    Eew = integrate.trapezoid(trough.incident_energy_on_system(alignment='EW', Ib=Ib) ).round(2)
    Etow = integrate.trapezoid(tower.incident_energy_on_system(Ib)).round(2)
    print("Total incident energy of NS Solar Trough ={} [MWh]".format(Ens))
    print("Total incident energy of NS Solar Trough ={} [MWh]".format(Eew))

    if save_to_file:
        IMG_FOLDER  = pathlib.Path("imgs/")
        IMG_FOLDER.mkdir(parents=True, exist_ok=True)
        plt.savefig(IMG_FOLDER/'plot-{}.png'.format(date.strftime('%m-%d')))
    # tower = solarII(df_irr.dni,1,IAM_tow(hoy),225000,99.3)
    # trough = di_sst(df_irr.dni,costhetai_NS(),IAM_tro(hoy),Tr, 5.76, 0.07, 18, 25, 1800)
# %%
