"""PFR isothermic energy balance module."""
from IPython.display import display

import numpy as np
from numpy.typing import NDArray

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Isothermic:
    r"""PFR isothermic energy balance.

    .. math::
        \frac{dT}{dz}=0

    .. math::
        \frac{dT_r}{dz}=0

    | :math:`T`: reactor's temperature.
    | :math:`T_r`: refrigerant's temperature.
    | :math:`z`: reactor's length coordinate.

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

    def initial_profile(self, reactor: PFR) -> NDArray[np.float64]:
        """Set initial energy profile in isothermic PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        NDArray[np.float64]
            initial temperature profile in all grid
        """
        return np.full(reactor.grid_size, self.temperature)

    def update_profile(
        self, reactor: PFR, variables: NDArray[np.float64]
    ) -> None:
        """Update PFR's temperature profile.

        Parameters
        ----------
        reactor : PFR
            PFR object.
        variables : NDArray[np.float64]
            Variables of solve_bvp ode solver.
        """
        reactor.temperature_profile = variables[-2, :]
        reactor.refrigerant_temperature_profile = None

    def border_conditions(self, reactor: PFR) -> tuple:
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

    def evaluate_balance(self, reactor: PFR) -> NDArray[np.float64]:
        """Evaluate energy balance.

        Parameters
        ----------
        reactor : PFR
          PFR object

        Returns
        -------
        NDArray[np.float64]
            Temperature rate of each substance on each reactor's z
            This is zero in a isothermic balance.
        """
        return np.zeros(reactor.grid_size)

    def __repr__(self) -> str:
        """Print LaTex representation."""
        latex1 = r"\frac{dT}{dz}=0"
        latex2 = r"\frac{dT_r}{dz}=0"
        return latex1, latex2
