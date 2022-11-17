import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealSolution(AbstractMix):
    def __init__(self, substance_list: list[Substance]):

        self.substances = substance_list

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        molar_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )

        total_molar_vol = np.dot(zi, molar_volumes)
        concentrations = np.divide(zi, total_molar_vol)
        return concentrations

    def volume(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )
        return np.dot(pure_volumes, zi)

    def mix_heat_capacity(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_cp = np.array(
            [
                substance.heat_capacity_liquid(temperature)
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
            enthalpies = np.append(enthalpies, substance.formation_enthalpy)

        return enthalpies

    def formation_enthalpies_correction(self, temperature: float, *args):

        enthalpies = np.array([])
        for substance in self.substances:
            if substance.normal_melting_point > 298.15:
                dhs = substance.heat_capacity_solid_dt_integral(
                    298.15, substance.normal_melting_point
                )
                dhf = substance.fusion_enthalpy(substance.normal_melting_point)
                dhl = substance.heat_capacity_liquid_dt_integral(
                    substance.normal_melting_point, temperature
                )

                enthalpies = np.append(enthalpies, dhs + dhf + dhl)
            else:
                enthalpies = np.append(
                    enthalpies,
                    substance.heat_capacity_liquid_dt_integral(
                        298.15, temperature
                    ),
                )
        return enthalpies
