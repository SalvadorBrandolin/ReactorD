from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isothermic:
    def __init__(self, temperature: float) -> None:
        self.temperature = temperature

    @property
    def irepr(self):
        display(symbols(self.__repr__()[0]))
        display(symbols(self.__repr__()[1]))

    def initial_profile(self, reactor: PFR):
        return np.full(reactor.grid_size, self.temperature)

    def update_profile(self, reactor: PFR, variables):
        reactor.temperature_profile = variables[-2, :]
        reactor.refrigerant_temperature_profile = None

    def border_conditions(self, reactor: PFR):
        return self.temperature, None

    def evaluate_balance(self, reactor: PFR):
        return np.zeros(reactor.grid_size)

    def __repr__(self) -> str:
        latex = (r"\frac{dT}{dz}=0", r"\frac{dT_r}{dz}=0")
        return latex
