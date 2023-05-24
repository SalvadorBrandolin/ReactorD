"""PFR adiabatic energy balance module."""
from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Adiabatic:
    """PFR Adiabatic energy balance class.

    Parameters
    ----------
    temperature_in_or_out : dict
    temperature of each substance in the inlet or outlet position.

    Raises
    ------
    NotImplementedError
    Both in and out border condition for temperature not implemented.
    """

    def __init__(self, temperature_in_or_out: dict) -> None:
        self.t_in_or_out = temperature_in_or_out

        if self.t_in_or_out.get("in") and self.t_in_or_out.get("out"):
            raise NotImplementedError(
                "Both in and out border condition for temperature not"
                " implemented."
            )

        if self.t_in_or_out.get("in"):
            self.t_value = self.t_in_or_out.get("in")
        else:
            self.t_value = self.t_in_or_out.get("out")

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()[0]))
        display(symbols(self.__repr__()[1]))

    def initial_profile(self, reactor: PFR):
        """Set initial energy profile in adiabatic PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            initial temperature profile in all grid
        """
        reactor.kinetic.set_dh_function()
        return np.full(reactor.grid_size, self.t_value)

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
        ndarray
            array with temperature in the inlet and outlet reactor
        """
        return self.t_in_or_out.get("in"), self.t_in_or_out.get("out")

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
        """
        delta_hs = reactor.kinetic.dhs_evaluate(
            reactor.temperature_profile, reactor.pressure_profile
        )
        cps = reactor.mix.mix_heat_capacity(
            reactor.mole_fraction_profile,
            reactor.temperature_profile,
            reactor.pressure_profile,
        )
        total_mole_flows = reactor.mass_profile.sum(axis=0)

        numerator = -np.multiply(delta_hs, reactor.r_rates_profile).sum(axis=0)
        denominator = np.multiply(cps, total_mole_flows)

        return np.divide(numerator, denominator) * reactor.transversal_area

    def __repr__(self) -> str:
        """Print LaTex representation."""
        latex1 = (
            r"\frac{dT}{dz}=\frac{\sum\Delta{H_j}_{(TP)}r_{j}}{{c_p}_{mi"
            r"x}{\sum}F_i}"
        )
        latex2 = r"\frac{dT_r}{dz}=0"
        return latex1, latex2
