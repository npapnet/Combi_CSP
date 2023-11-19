Using CombiCSP for Solar Power Calculations in Crete
====================================================

This tutorial demonstrates how to use the CombiCSP package to calculate solar power for a project in Crete. We will use both Solar Tower and Solar Trough models.

Prerequisites
-------------

Before starting, ensure you have the following libraries installed:

- matplotlib
- scipy
- pandas
- ipykernel
- pvlib_python
- iapws
- numpy-financial

Import Necessary Libraries
--------------------------

First, import the required modules:

.. code-block:: python

    import pathlib
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.interpolate import make_interp_spline, BSpline
    from CombiCSP.misc import heatmap2d
    from CombiCSP import SolarTroughCalcs, SolarTowerCalcs, HOYS_DEFAULT, SolarSystemLocation
    from CombiCSP.storage import Tr

Load Data and Constants
-----------------------

Load the hourly data and define the location for the solar calculations:

.. code-block:: python

    hoy = HOYS_DEFAULT
    FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
    df_pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'])
    Ib = df_pvgis.loc[:,'Gb(n)']
    sslCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)

Calculate Tower Power Output
----------------------------

Set up the tower dimensions and calculate the power output:

.. code-block:: python

    stc = SolarTowerCalcs(alt = 200*10e-3, Ht = 0.1, Ar = 99.3, A_helio = 225000, slobj=sslCrete)
    oTow = stc.perform_calc(Ib)

Plot Tower Output
-----------------

Plot the power output from the tower:

.. code-block:: python

    plt.plot(hoy, oTow.data, label='1')
    plt.xlabel('Time (hour of year)')
    plt.ylabel('Power (MW)')
    plt.title('Tower')
    plt.legend()

Calculate Trough Power Output
-----------------------------

Set up the trough dimensions and calculate the power output for both East-West and North-South orientations:

.. code-block:: python

    sotr = SolarTroughCalcs(foc_len = 0.88, N = 1800, L = 25, Ws = 18, Wr = 0.07, Wc = 5.76, slobj=sslCrete)
    oew = sotr.perform_calcs_EW(Ib=Ib, Tr=Tr)
    ons = sotr.perform_calcs_NS(Ib=Ib, Tr=Tr)

Plot Trough Output
------------------

Visualize the power output from the troughs:

.. code-block:: python

    plt.plot(hoy, ons.data)
    plt.plot(hoy, oew.data)
    plt.xlabel('Time (hour of year)')
    plt.ylabel('Power (MW)')
    plt.legend(('EW','NS'))
    plt.title('Trough')
    plt.show()

Generate Heatmaps
-----------------

Finally, generate heatmaps for the tower and troughs:

.. code-block:: python

    plt.title('Tower')
    heatmap2d(oTow.data4surf().T)
    plt.title('Trough N-S')
    heatmap2d(ons.data4surf().T)
    plt.title('Trough E-W')
    heatmap2d(oew.data4surf().T)
