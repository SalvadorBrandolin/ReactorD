"""Ideal gas Module."""
from typing import List

import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealGas(AbstractMix):
    """IdealGas object class.

    This class should be used when the mixture is considered an ideal gas.

    Parameters
    ----------
    substance_list : list [float]
        list of substance objects

    Attributes
    ----------
    substances : list [float]
        list of substance objects
    """

    def __init__(self, substance_list: List[Substance]):
        self.substances = substance_list

    def concentrations(self, moles, temperature, pressure):
        """Calculate concentrations of the mixture.

        Parameters
        ----------
        moles : ndarray or list [float]
            Moles of each substance
        temperature : float
            System temperature
        pressure : float
            System pressure

        Returns
        -------
        ndarray or list [float]
            Concentration of each substance
        """
        zi = self.mol_fracations(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        density = pressure / (r * temperature)
        return np.multiply(zi, density)

    def volume(self, moles, temperature, pressure):
        """Calculate the volume of the mixture.

        Parameters
        ----------
        moles : ndarray or list[float]
            Moles of each substance
        temperature : float
            System temperature
        pressure : float
            System pressure

        Returns
        -------
        float
            Volume of the mixture
        """
        total_moles = np.sum(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = total_moles * r * temperature / pressure
        return volume

    def mix_heat_capacity(self, moles, temperature, *args):
        """Calculate heat capacity of th mixture.

        Parameters
        ----------
        moles : ndarray or list [float]
            Moles of each substance
        temperature : float
            System temperature

        Returns
        -------
        float
            Heat capacity of the mixture
        """
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
        """Return the ideal gas formation enthalpies in a ordered ndarray.

        Method that read the ideal gas formation enthalpies of mix
        and returns them in a ordered ndarray.

        Returns
        -------
        ndarray [float]
            Ideal gas formation enthalpies of each substance
        """
        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(enthalpies, substance.formation_enthalpy_ig)

        return enthalpies

    def formation_enthalpies_correction(self, temperature: float, *args):
        """Calculate the enthalpy of formation at the specified temperature.

        Method that corrects the enthalpy of formation of pure
        substances of 298.15 K 100000 Pa at the specified temperature
        using Kirchhoff's equation.

        Parameters
        ----------
        temperature : float
            Correction temperature for the formation enthalpies. [K]
        """
        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(
                enthalpies,
                substance.heat_capacity_gas_dt_integral(298.15, temperature),
            )

        return enthalpies
