import numpy as np
from scipy.integrate import quad
from abc import ABCMeta 
from Mix import Abstract_Mix
from kinetics import Kinetics

class ReactorBase(metaclass=ABCMeta):

    def __init__(
        self, 
        mix : Abstract_Mix,
        list_of_reactions : dict[str:function],
        stoichiometry : list,
        **options,
    ):                      
        
        self.kinetic : Kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            **options)
        self.stoichiometry = np.array(stoichiometry)
        self.mix = mix
        self.list_of_reactions = list_of_reactions
        
        
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