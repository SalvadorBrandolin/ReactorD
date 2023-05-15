from typing import Callable

import numpy as np
from numpy.typing import NDArray

from reactord.mix.abstract_mix import AbstractMix
from reactord.flowreactors.stationary_1d.pfr.pfr import PFR
from .noisothermic_all_constant import NoIsothermicAllConstant


class NoIsothermicWithUConstant(NoIsothermicAllConstant):
    def __init__(
        self,
        temperature_in_or_out: dict,
        refrigerant_mix: AbstractMix,
        refrigerant_in_temperature: float,
        heat_exchange_coefficient: Callable,
        exchange_operation: str = "cocurrent",
    ) -> None:
        super().__init__(
            temperature_in_or_out,
            refrigerant_in_temperature,
            heat_exchange_coefficient,
            exchange_operation,
        )

        self.ref_mix = refrigerant_mix

    def evaluate_balance(self, reactor: PFR) -> NDArray:
        # =====================================================================
        # Reactor energy balance
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

        numerator = self.operation_coeff * self.u * a * (
            reactor.refrigerant_temperature_profile
            - reactor.temperature_profile
        ) - np.multiply(delta_hs, reactor.r_rates_profile).sum(axis=0)

        denominator = np.multiply(cps, total_mole_flows)

        dt_dz = np.divide(numerator, denominator)

        # =====================================================================
        # Refrigerant energy balance
        # =====================================================================
        dta_dz = np.zeros(reactor.grid_size)

        return np.vstack(dt_dz, dta_dz)
