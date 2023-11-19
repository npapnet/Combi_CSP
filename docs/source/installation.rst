Installation/Setup
==================

Requirements
------------

The package depends on the following "standard" libraries:

- matplotlib
- scipy
- pandas
- ipykernel

Additional libraries required are:

- pvlib_python
- iapws (The International Association for the Properties of Water and Steam)
- numpy-financial

Conda Installation
------------------

The following describes a minimum environment setup using Conda.

(Optional) It is preferable to create a new environment for the packages:

.. code-block:: bash

    conda create -n combicsp python=3

The installation requires:

- matplotlib
- scipy
- pandas
- ipykernel

To install these packages, use the command:

.. code-block:: bash

    conda install matplotlib scipy pandas ipykernel

For the additional libraries:

- pvlib_python:

  .. code-block:: bash

      conda install -c conda-forge pvlib-python

- iapws (The International Association for the Properties of Water and Steam):

  .. code-block:: bash

      conda install -c conda-forge iapws

- numpy-financial:

  .. code-block:: bash

      conda install -c conda-forge numpy-financial
