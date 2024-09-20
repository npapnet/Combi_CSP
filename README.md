# Combi_CSP
CombiCSP is an open source software for dynamic modelling of concentrating solar energy power plants. CombiCSP utilizes solar resource, system engineering inputs as well as financial tools to provide dynamic simulations and annual yields of concentrating solar power plants. It readily provides modelling of plants based on solar power tower and parabolic trough collectors and it can be extended to novel solar energy modeling approaches and analyses as needed.

The relevant work pertaining to this software is : https://doi.org/10.1016/j.apenergy.2022.119450

DOCUMENTATION: [https://npapnet.github.io/Combi_CSP/](https://npapnet.github.io/Combi_CSP/)

# Quick-Start

There are several methods to install 
## Using pip

in a new environmetn 

`pip  install CombiCSP`

## Cloning the repository and using setuptools

clone the repository from [the CombiCSP][2].

then move to the root of the repository and perform (this requires python [`setuptools`][1] )

`python setup.py install`

The setup should take care of all problems.

The `Combi_CSP_oop.ipynb` describes a typical use scenario.

Additionally, example cases are scripted in the following files:

- `CSPCret.py`: heliostat and CSP power and energy outputs in a location in Crete, Greece
- `CSP50Compare.py`: Combined heliostat and CSP power and energy outputs in a location in Crete, Greece

##  Cloning the repository and Manual installation/setup

### requirements

The following **mainstream** packages are required for this library (most of them are already installed in a typical installation).

- matplotlib
- scipy
- pandas
- ipykernel
- seaborn

Additional libraries are:

- pvlib_python
- iapws (The InternationalAssociation for the Properties of Water and Steam)
- numpy-financial

### Conda installation

The following describes a minimum environment using conda. 

(Optional) Preferably create a new environment for the packages

`conda create -n combicsp python=3`

installation requires:

- matplotlib
- scipy
- pandas
- ipykernel

`conda install matplotlib scipy pandas ipykernel`

- pvlib_python:

`conda install -c conda-forge pvlib-python`

- iapws: The InternationalAssociation for the Properties of Water and Steam

`conda install -c conda-forge iapws`

- numpy-financial

`conda install -c conda-forge numpy-financial`

# Developers

- G. E. Arnaoutakis: Technical direction, and original calculation functions
- N. Papadakis: Software Design, Package Maintenance.


# TODO

See [history.md](history.md)

[1]: https://pypi.org/project/setuptools/
[2]: https://github.com/npapnet/Combi_CSP.git