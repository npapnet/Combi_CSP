
# Scope

The scope of this file is to describe the contents of the examples folder.

## `CSPCrete.py` 

TODO provide description

## `CSP50Compare_combi.py` 

TODO Provide description

## `CSP50Compare_np.py` 

TODO Provide description

## `CSP50common_econ.py`

Contains investment parameters for the economic analysis of the Combined Solar Power 50MW  project.

TODO Provide description

## `Tower_modelling.py`

This is a file created on the 18th of September 2024 to model only the tower of the CSP plant. It is a work in progress.

It is primarily done to remember what I have done so far, and see what I can improve in terms of the code.

Things I have idntified so far:

- Better to move to a dataframe for the time series. 
 
  -  I have extended the output container, in a way that does not affect the test, but I need to make this more efficient. 
  -  I should remove functions like `hour_power_arr()` in favour of the dataframe

-  I have created a class diagram.