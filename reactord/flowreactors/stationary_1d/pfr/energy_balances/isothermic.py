import numpy as np

from IPython.display import display

from sympy import symbols

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR


class Isothermic:
    def __init__(self, temperature: float) -> None:
        self.temperature = temperature
        
    @property
    def irepr(self):
        display(symbols(self.__repr__()[0]))
        display(symbols(self.__repr__()[1]))

    def initial_profile(self, reactor: PFR):
        temperature_profile = np.full(reactor.grid_size, self.temperature)
        refrigerant_temperature_profile = np.zeros(reactor.grid_size)
        return temperature_profile, refrigerant_temperature_profile

    def border_conditions(self, reactor: PFR):
        return np.array([self.temperature, 0]), np.array([None, None])

    def evaluate_balance(self, reactor: PFR):
        return np.zeros(reactor.grid_size), np.zeros(reactor.grid_size)

    def __repr__(self) -> str:
        latex = (r"\frac{dT}{dz}=0", r"\frac{dT_r}{dz}=0")
        return latex
