import numpy as np
from reactord.substance import Substance
from .abstract_mix import AbstractMix


class IdealSolution(AbstractMix):
    def __init__(self, substance_list: list[Substance]):

        self.substances = substance_list
        self.formation_enthalpies = self.formation_enthalpies_correction()

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        molar_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )

        total_molar_vol = np.dot(zi, molar_volumes)
        concentrations = np.divide(zi, total_molar_vol)  # moles/m^3
        return concentrations

    def volume(self, moles, temperature, pressure):
        pure_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )
        return np.dot(pure_volumes, moles)

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

    def formation_enthalpies_correction(self):
        # for substance in self.substances:
        enthalpies = np.array([])
        for substance in self.substances:
            if substance.normal_melting_point > 298.15:
                dhs = substance.heat_capacity_solid_dt_integral(
                    298.15, substance.normal_melting_point
                )
                dhf = substance.fusion_enthalpy(substance.normal_melting_point)
                dhl = substance.heat_capacity_liquid_dt_integral(
                    substance.normal_melting_point, 298.15
                )
                print(dhs)
                print(dhf)
                print(dhl)
                print(substance.formation_enthalpy)

                enthalpies = np.append(
                    enthalpies, substance.formation_enthalpy + dhs + dhf + dhl
                )
            else:
                enthalpies = np.append(enthalpies, substance.formation_enthalpy)
        return enthalpies
