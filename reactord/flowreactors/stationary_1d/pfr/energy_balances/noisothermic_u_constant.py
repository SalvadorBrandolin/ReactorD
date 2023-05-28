"""PFR noisothermic energy balance module."""
import numpy as np
from numpy.typing import NDArray

from reactord.flowreactors.stationary_1d.pfr.energy_balances import (
    NoIsothermicAllConstant,
)
from reactord.flowreactors.stationary_1d.pfr.pfr import PFR
from reactord.mix.abstract_mix import AbstractMix


class NoIsothermicUConstant(NoIsothermicAllConstant):
    r"""PFR noisothermic energy balance.

    Heat exchange coefficient U remains constant. Refrigerant's pressure
    profile, composition profile and volumetric flow are assumed constant.

    .. math::
        \frac{1}{a_t}\frac{dT}{dz}=\frac{Ua(T_r-T)+\sum\Delta{H_j}r_{j}}{{c_p}_
        {mix}{\sum}F_i}

    .. math::
        \text{Co-current refrigerant energy balance: }
        \frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T-T_r)}{{c_p}_{r}{\sum}F_i}

    .. math::
        \text{Counter-current refrigerant energy balance: }
        \frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T_r-T)}{{c_p}_{r}{\sum}F_i}

    | :math:`a_t`: reactor's transversal area.
    | :math:`U`: heat exchange coefficient.
    | :math:`a`: reactor's heat exchange area per reactor's volume.
    | :math:`T`: reactor's temperature.
    | :math:`T_r`: refrigerant's temperature.
    | :math:`z`: reactor's length coordinate.
    | :math:`\Delta{H_j}_{(TP)}`: enthalpy of `j`-th reaction.
    | :math:`r_{j}`: rate of `j`-th reaction.
    | :math:`{c_p}_{mix}`: mixture's heat capacity.
    | :math:`{c_p}_{r}`: refrigerant's heat capacity.
    | :math:`F_i`: molar flow of the `i`-th substance.

    Parameters
    ----------
    temperature_in_or_out: dict
        Temperature of reactive mixture in the inlet or outlet position. Use
        {"in": 200} for 200 K at inlet or {"out": 200} for 200 K at outlet.
    refrigerant_in_temperature: float
        Refrigerant temperature. [K]
    heat_exchange_coefficient: float
        Heat exchange coefficient. [J/K/s/m²]
    refrigerant_mix: AbstractMix
        Refrigerant mixture.
    efrigerant_composition: NDArray[np.float64]
        Mole fractions of the refrigerant mixture substances.
    refrigerant_pressure: float
        Refrigerant's pressure. [Pa]
    refrigerant_moalr_flow: float
        Refrigerant's molar flow. [mol/s]
    heat_exchange_arrangement: str, optional
        Options: 'co-current', 'counter-current', by default 'co-current'.
    """

    def __init__(
        self,
        temperature_in_or_out: dict,
        refrigerant_in_temperature: float,
        heat_exchange_coefficient: float,
        refrigerant_mix: AbstractMix,
        refrigerant_composition: NDArray[np.float64],
        refrigerant_pressure: float,
        refrigerant_molar_flow: float,
        heat_exchange_arrangement: str = "co-current",
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
        self.refrigerant = refrigerant_mix
        self.refrigerant_composition = refrigerant_composition
        self.refrigerant_pressure = refrigerant_pressure
        self.refrigerant_flow = refrigerant_molar_flow
        self.arrangement = heat_exchange_arrangement

        if self.arrangement == "co-current":
            self.arrangement_coeff = 1
        elif self.arrangement == "counter-current":
            self.arrangement_coeff = -1
        else:
            raise ValueError(
                f"{self.arrangement} is not a valid heat exchange"
                "arrangement."
            )

    def border_conditions(self, reactor: PFR) -> tuple:
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        temperature at in and out condition.
        """
        if self.arrangement == "co-current":
            in_border_cond = np.array(
                [self.t_in_or_out.get("in"), self.refrigerant_in_temperature]
            )
            out_border_cond = np.array([self.t_in_or_out.get("out"), None])
        else:
            in_border_cond = np.array([self.t_in_or_out.get("in"), None])
            out_border_cond = np.array(
                [self.t_in_or_out.get("out"), self.refrigerant_in_temperature]
            )

        return in_border_cond, out_border_cond

    def evaluate_balance(self, reactor: PFR) -> NDArray[np.float64]:
        """Evaluate energy balance.

        Parameters
        ----------
        reactor : PFR
          PFR object

        Returns
        -------
        NDArray[np.float64]
            Temperature rate of each substance on each reactor's z.
        """
        # =====================================================================
        # reactor's energy balance
        # =====================================================================
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

        # =====================================================================
        # refrigerant's energy balance
        # =====================================================================
        ref_pressure = np.full(reactor.grid_size, self.refrigerant_pressure)

        cps_ref = self.refrigerant.mix_heat_capacity(
            self.refrigerant_composition,
            reactor.refrigerant_temperature_profile,
            ref_pressure,
        )

        numerator_ref = (
            self.arrangement_coeff
            * self.u
            * a
            * (
                reactor.temperature_profile
                - reactor.refrigerant_temperature_profile
            )
        )

        denominator_ref = self.refrigerant_flow * cps_ref

        dta_dz = (
            np.divide(numerator_ref, denominator_ref)
            * reactor.transversal_area
        )

        return np.vstack((dt_dz, dta_dz))

    def __repr__(self) -> str:
        """Print LaTex representation."""
        latex1 = (
            r"\frac{1}{a_t}\frac{dT}{dz}=\frac{Ua(T_a-T)+\sum\Delta{H_j}r_{j}}"
            r"{{c_p}_{mix}{\sum}F_i}"
        )

        if self.arrangement == "co-current":
            latex2 = (
                r"\frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T-T_r)}{{c_p}_{r}{\sum"
                r"}F_i}}"
            )
        else:
            latex2 = (
                r"\frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T_r-T)}{{c_p}_{r}{\sum"
                r"}F_i}}"
            )

        return latex1, latex2
