from typing import List

import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealGas(AbstractMix):
    def __init__(self, substance_list: List[Substance]):
        self.substances = substance_list

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        density = pressure / (r * temperature)
        return np.multiply(zi, density)

    def volume(self, moles, temperature, pressure):
        total_moles = np.sum(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = total_moles * r * temperature / pressure
        return volume

    def mix_heat_capacity(self, moles, temperature, *args):
        zi = self.mol_fracations(moles)
        pure_cp = np.array(
            [
                substance.heat_capacity_gas(temperature)
                for substance in self.substances
            ]
        )
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp

    def _formation_enthalpies_set(self):
        """Method that read the ideal gas formation enthalpies of mix's
        and returns them in a ordered ndarray.
        """
        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(enthalpies, substance.formation_enthalpy_ig)

        return enthalpies

    def formation_enthalpies_correction(self, temperature: float, *args):
        """Method that correction therm for the formation enthalpy of
        the pure substances from 298.15 K 100000 Pa to temperature and
        pressure using the eq:

        Parameters
        ----------
        temperature : float
            Correction temperature for the formation enthalpies. [K]
        pressure : float
            Correction pressure for the formation enthalpies. [Pa]

        """

        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(
                enthalpies,
                substance.heat_capacity_gas_dt_integral(298.15, temperature),
            )

        return enthalpies
