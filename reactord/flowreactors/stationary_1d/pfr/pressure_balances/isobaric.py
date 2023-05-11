from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isobaric:
    def __init__(self, pressure: float) -> None:
        self.pressure = pressure

    @property
    def irepr(self):
        display(symbols(self.__repr__()))

    def initial_profile(self, reactor: PFR):
        return np.full(reactor.grid_size, self.pressure)

    def update_profile(self, reactor: PFR, variables):
        reactor.pressure_profile = variables[-1, :]

    def border_conditions(self, reactor: PFR):
        return self.pressure, None

    def evaluate_balance(self, reactor: PFR):
        return np.zeros(reactor.grid_size)

    def __repr__(self):
        return r"\frac{dP}{dz}=0"
