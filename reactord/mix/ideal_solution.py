"""Ideal solution Module."""
from typing import List

import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealSolution(AbstractMix):
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
        moles : ndarray or list[float]
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
        """Calculate the volume of the mixture.

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
        float
            Volume of the mixture
        """
        zi = self.mol_fracations(moles)
        pure_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )
        return np.dot(pure_volumes, zi)

    def mix_heat_capacity(self, moles, temperature, pressure):
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
                substance.heat_capacity_liquid(temperature)
                for substance in self.substances
            ]
        )
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp

    def _formation_enthalpies_set(self):
        """Return the formation enthalpies in a ordered ndarray.

        Method that read the formation enthalpies of mix and returns
        them in a ordered ndarray.

        Returns
        -------
        ndarray [float]
            Formation enthalpies of each substance
        """
        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(enthalpies, substance.formation_enthalpy)

        return enthalpies

    def formation_enthalpies_correction(self, temperature: float, *args):
        """Return corrects the enthalpy of formation of pure substances.

        Method that corrects the enthalpy of formation of pure
        substances When its melting temperature is greater than 298.

        """
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
