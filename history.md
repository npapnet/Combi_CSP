# Revision History/Versions

# TODO items

- HOYS_DEFAULT should be removed. The Ib data has a time series and that should be used instead of HOYS.

- Add description to the SolarTrough and SolarTower
- Rename calculation functions solarII (Tower) and di_sst (Trough), and add all the relebant parameters
- Remove/rewrite the following files:
  - CSP.py
  - storage.py
  - Transmittance.py
- (low priority) Documentation
- (low priority) submit to pypi

## Work in progress of v1.6.0 
### **(Working on)**
- Develop a new function for the combined system (CombiCSP_System).
- Develop unit tests for the combined system CombiCSP_System.
### **Completed**
- complete transition of `CSP50Compare.py`: there are financial optimisation function that need completion.
- `CSP50_common_econ.py` contains important information.
- Wrote an equivalent version of CSP50Compare.py using the new class system

- removed old procedural function of solar geometry
- removed old procedural functions that are now part of SolarTowerCalcs and SolarTroughCalcs
- removed original versions of:
  - `CSPCret.py` and comparison file
  - `CSP50Compare.py` and comparison file
  - `CombiCSP.ipynb`

- Removed the following files:
  - SolarGeometry.py


## v1.5.0 (20220721)

- created classes:
  - `SolarTowerCalcs`: Solar tower calculations 
  - `SolarTroughCalcs`: Solar trough calculations with EW and NS orientations
  - `SolarSystemLocation`:  performs various calculations regarding Solar Geometry at a specific location
- Now the basic calculations can be performed with classes in a more consise manner.
- It is now possible to set the location of the CSP installation (in the original version Crete/Ierapetra was the default and only installation site).
-  moved 
   - `CSPCret.py` to examples and created OOP equivalent and comparison
   - `CSP50Compare.py` to examples and created OOP equivalent and comparison
   - `CombiCSP.ipynb` to examples
- Create a  `CombiCSP-Quickstart.ipynb` to showcase the package functionality

## v1.0.0

This refers to the initial version of the package which can be installed with setup.py.