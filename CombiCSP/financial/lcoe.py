#%%
import numpy as  np
import pandas as pd

class LCOECalculator():
    def __init__(self):
        pass
        
    def lcoe_Eu_MWh(self, df:pd.DataFrame, discount_rate:float):
        self.discount_rate = discount_rate
        # Calculate the discounted costs over the lifetime of the project
        self.discounted_costs = df['costs_per_annum'] / (1 + self.discount_rate) ** df['year']
        # Calculate the total energy produced over the lifetime of the project
        self.discounted_energy = df['energy_MWh'] / (1 + self.discount_rate) ** df['year']
        # Calculate the LCOE
        lcoe = sum(self.discounted_costs) / sum(self.discounted_energy)
        return lcoe
# can you write a test for this function?

#%%
if __name__ == "__main__":
    year  =  range(0, 31)
    energy_MWh = np.ones(31) * (1600)
    costs_per_annum = np.ones(31) * 10000
    costs_per_annum[0] = 1e6

    df = pd.DataFrame({'year': year,
        'energy_MWh':energy_MWh, 
        'costs_per_annum': costs_per_annum
        })    
        
    # %%
    discount_rate = 0.04
    discounted_costs = df['costs_per_annum'] / (1 + discount_rate) ** df['year']
        
    # %%

    1 / (1 + discount_rate) ** df['year']
        
    # %%
    lcoecalc = LCOECalculator() 
    lcoecalc.lcoe_Eu_MWh(df, 0.1)
    # %%
