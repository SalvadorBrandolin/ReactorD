"""PFR noisothermic energy balance module."""
from IPython.display import display

import numpy as np
from numpy.typing import NDArray

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class NoIsothermicAllConstant:
    """PFR noisothermic energy balance class.

    Parameters
    ----------
    temperature_in_or_out: dict
    refrigerant_in_temperature: float
    heat_exchange_coefficient: float
    """

    def __init__(
        self,
        temperature_in_or_out: dict,
        refrigerant_in_temperature: float,
        heat_exchange_coefficient: float,
    ) -> None:
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

        self.refrigerant_in_temperature = refrigerant_in_temperature
        self.u = heat_exchange_coefficient

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()[0]))
        display(symbols(self.__repr__()[1]))

    def initial_profile(self, reactor: PFR) -> NDArray:
        """Set initial energy profile in noisothermic PFR.

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
        initial_profile = np.array(
            [
                np.full(reactor.grid_size, self.t_value),
                np.full(reactor.grid_size, self.refrigerant_in_temperature),
            ]
        )
        return initial_profile

    def update_profile(self, reactor: PFR, variables: NDArray) -> None:
        """Update profile."""
        reactor.temperature_profile = variables[-3, :]
        reactor.refrigerant_temperature_profile = variables[-2, :]

    def border_conditions(self, reactor: PFR):
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        temperature at in and out condition
        """
        in_border_cond = np.array(
            [self.t_in_or_out.get("in"), self.refrigerant_in_temperature]
        )
        out_border_cond = np.array([self.t_in_or_out.get("out"), None])

        return in_border_cond, out_border_cond

    def evaluate_balance(self, reactor: PFR) -> NDArray:
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

        a = 4 / (2 * reactor.tube_radius)

        numerator = self.u * a * (
            reactor.refrigerant_temperature_profile
            - reactor.temperature_profile
        ) - np.multiply(delta_hs, reactor.r_rates_profile).sum(axis=0)
        denominator = np.multiply(cps, total_mole_flows)

        dt_dz = np.divide(numerator, denominator) * reactor.transversal_area
        dta_dz = np.zeros(reactor.grid_size)

        return np.vstack((dt_dz, dta_dz))

    def __repr__(self) -> str:
        """Print LaTex representation."""
        latex1 = (
            r"\frac{dT}{dz}=\frac{Ua(T_a-T)+\sum\Delta{H_j}r_{j}}{{c_p}_{mix}"
            r"{\sum}F_i}"
        )
        latex2 = r"\frac{dT_r}{dz}=0"

        return latex1, latex2
