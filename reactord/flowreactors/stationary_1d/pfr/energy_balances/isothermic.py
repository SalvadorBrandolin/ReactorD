"""PFR isothermic energy balance module."""
from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isothermic:
    """PFR isothermic energy balance class.

    Parameters
    ----------
    temperature: float
        Isothermic temperature
    """

    def __init__(self, temperature: float) -> None:
        self.temperature = temperature

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()[0]))
        display(symbols(self.__repr__()[1]))

    def initial_profile(self, reactor: PFR):
        """Set initial energy profile in isothermic PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            initial temperature profile in all grid
        """
        return np.full(reactor.grid_size, self.temperature)

    def update_profile(self, reactor: PFR, variables):
        """Update profile."""
        reactor.temperature_profile = variables[-2, :]
        reactor.refrigerant_temperature_profile = None

    def border_conditions(self, reactor: PFR):
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        isothermic temperature
        """
        return self.temperature, None

    def evaluate_balance(self, reactor: PFR):
        """Evaluate energy balance.

        Parameters
        ----------
        reactor : PFR
          PFR object

        Returns
        -------
        ndarray
            Temperature rate of each substance on each reactor's z
            This is zero in a isothermic balance.
        """
        return np.zeros(reactor.grid_size)

    def __repr__(self) -> str:
        latex1 = r"\frac{dT}{dz}=0"
        latex2 = r"\frac{dT_r}{dz}=0"
        return latex1, latex2
