"""PFR isobaric pressure balance module."""
from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isobaric:
    """PFR isobaric energy pressure class.

    Parameters
    ----------
    pressure : float
        Isobaric pressure
    """

    def __init__(self, pressure: float) -> None:
        self.pressure = pressure

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()))

    def initial_profile(self, reactor: PFR):
        """Set initial pressure profile in non-isobaric PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            initial pressure profile in all grid
        """
        return np.full(reactor.grid_size, self.pressure)

    def update_profile(self, reactor: PFR, variables):
        """Update profile."""
        reactor.pressure_profile = variables[-1, :]

    def border_conditions(self, reactor: PFR):
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            array with pressure in the inlet and outlet reactor
        """
        return self.pressure, None

    def evaluate_balance(self, reactor: PFR):
        """Evaluate energy balance.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            Isobaric pressure of each substance on each reactor's z
        """
        return np.zeros(reactor.grid_size)

    def __repr__(self):
        """Represent equation of PFR isobaric pressure balance."""
        return r"\frac{dP}{dz}=0"
