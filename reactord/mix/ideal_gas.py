import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealGas(AbstractMix):
    def __init__(self, substance_list: list[Substance]):
        self.substances = substance_list
        self.formation_enthalpies = [
            substance.formation_enthalpy_ig for substance in self.substances
        ]

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)

        molar_volumes = np.array(
            [
                substance.volume_gas(temperature, pressure)
                for substance in self.substances
            ]
        )

        total_molar_vol = np.dot(zi, molar_volumes)
        concentrations = np.divide(zi, total_molar_vol)  # moles/m^3
        return concentrations

    def volume(self, moles, temperature, pressure):
        pure_volumes = np.array(
            [
                substance.volume_gas(temperature, pressure)
                for substance in self.substances
            ]
        )
        return np.dot(pure_volumes, moles)

    def mix_heat_capacity(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_cp = np.array(
            [
                substance.heat_capacity_gas(temperature)
                for substance in self.substances
            ]
        )
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp

    def partial_pressures(self, moles, temperature, pressure):
        """method that calculates the partial pressures of the mixture

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]

        Returns
        -------
        ndarray
            array that contains the partial pressures of mixture's
            substances
        """
        zi = self.mol_fracations(moles)
        partial_pressures = np.multiply(zi, pressure)
        return partial_pressures

    def partial_p_to_conc(self, partial_pressures, temperature):
        r = 8.31446261815324  # J/mol.K
        self.partial_pressures = np.array(partial_pressures)
        conc = self.partial_pressures / (r * temperature)  # mol/m^3
        return conc
