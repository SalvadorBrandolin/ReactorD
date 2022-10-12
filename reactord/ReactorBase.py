import numpy as np
from scipy.integrate import quad
from abc import ABCMeta, abstractmethod 
from Mix import Abstract_Mix
from kinetics import Kinetics

class ReactorBase(metaclass=ABCMeta):

    def __init__(
        self, 
        mix: Abstract_Mix,
        list_of_reactions: list[function],
        stoichiometry: list[float],
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
        pass

    @abstractmethod
    def _border_condition_builder(self):
        pass

    @abstractmethod
    def _mass_balance(self):
        pass

    @abstractmethod
    def _energy_balance(self):
        pass

    @abstractmethod
    def _pressure_balance(self):
        pass

    def _particle_mass_balance(self):
        pass
    
    def _particle_energy_balance(self):
        pass

    @abstractmethod
    def solve(self):
        pass