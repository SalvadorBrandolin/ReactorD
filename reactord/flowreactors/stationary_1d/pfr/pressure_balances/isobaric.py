"""PFR isobaric pressure balance module."""
from IPython.display import display

import numpy as np
from numpy.typing import NDArray

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isobaric:
    r"""PFR isobaric energy pressure class.

    .. math::
        \frac{dP}{dz}=0

    | :math:`P`: reactor's pressure.
    | :math:`z`: reactor's length coordinate.

    Parameters
    ----------
    pressure : float
        Isobaric pressure.
    """

    def __init__(self, pressure: float) -> None:
        self.pressure = pressure

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()))

    def initial_profile(self, reactor: PFR) -> NDArray[np.float64]:
        """Set initial pressure profile in non-isobaric PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object.

        Returns
        -------
        NDArray[np.float64]
            initial pressure profile in all grid.
        """
        return np.full(reactor.grid_size, self.pressure)

    def update_profile(
        self, reactor: PFR, variables: NDArray[np.float64]
    ) -> None:
        """Update th reactor's presure profile.

        Parameters
        ----------
        reactor : PFR
            PFR object.
        variables : NDArray[np.float64]
            Variables of solve_bvp ode solver.
        """
        reactor.pressure_profile = variables[-1, :]

    def border_conditions(self, reactor: PFR) -> tuple:
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        tuple
            array with pressure in the inlet and outlet reactor
        """
        return self.pressure, None

    def evaluate_balance(self, reactor: PFR) -> NDArray[np.float64]:
        """Evaluate pressure balance.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        NDArray[np.float64]
            Isobaric pressure of each substance on each reactor's z
        """
        return np.zeros(reactor.grid_size)

    def __repr__(self) -> str:
        """Represent equation of PFR isobaric pressure balance."""
        return r"\frac{dP}{dz}=0"
