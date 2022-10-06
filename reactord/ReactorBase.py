import numpy as np
from scipy.integrate import quad
from abc import ABCMeta 
from Mix import Abstract_Mix

class ReactorBase(metaclass=ABCMeta):

    def __init__(
        self, 
        mix : Abstract_Mix,
        dict_of_reactions : dict[str:function],
        stoichiometry : list,
        **options,
    ):                      
        
        if np.ndim(np.shape(stoichiometry)) == 1:
            self.num_reactions = 1
            self.total_substances = np.shape(stoichiometry)[0]
        else:
            self.num_reactions, self.total_substances = np.shape(stoichiometry) 

        if self.num_reactions != len(dict_of_reactions):
            raise IndexError(
                "'stoichiometry' rows number must be equal to" 
                "dict_of_reactions' length" 
            )
        
        if len(mix) != self.total_substances:
            raise IndexError(
                "'stoichiometry' columns number must be equal to substances" 
                "number in 'mix' object" 
            )

        self.stoichiometry = np.array(stoichiometry)
        self.mix = mix
        self.dict_of_reactions = dict_of_reactions
        self.std_reaction_enthalpies = np.dot(
            self.stoichiometry, self.mix.h_formations
        )
        
    def reaction_enthalpies(self, temperature, pressure):
        t_0 = 298.15
        pure_cp_integrals = self.mix.pure_heat_capacities_integrals(
            temperature, pressure
        )
        enthalpies_change = np.dot(self.stoichiometry, pure_cp_integrals)
        return (self.std_reaction_enthalpies + enthalpies_change)

    def _border_condition_builder(self):
        pass

    def _mass_balance(self):
        pass

    def _energy_balance(self):
        pass

    def _pressure_balance(self):
        pass

    def _particle_mass_balance(self):
        pass
    
    def _particle_energy_balance(self):
        pass
    
    def _result_report(self):
        pass

    def solve(self):
        pass