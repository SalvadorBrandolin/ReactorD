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

    def __init__(self, **substance_dict) -> None:

        substance_list = [
            value
            if isinstance(value, Substance)
            else Substance.from_thermo_database(value)
            for value in substance_dict.values()
        ]
        self.substances = substance_list

    def concentrations(
        self, moles: List[float], temperature: float, pressure: float
    ) -> List[float]:
        """Calculate concentrations of the mixture.

        Parameters
        ----------
        moles : ndarray or list[float]
            Moles of each substance
        temperature : float
            System temperature [T]
        pressure : float
            System pressure [Pa]

        Returns
        -------
        ndarray or list [float]
            Concentration of each substance [mol/m³]
        """
        mol_fractions = self.mol_fractions(moles)

        molar_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )

        total_molar_vol = np.dot(mol_fractions, molar_volumes)
        concentrations = np.divide(mol_fractions, total_molar_vol)

        return concentrations

    def volume(
        self, moles: List[float], temperature: float, pressure: float
    ) -> float:
        """Calculate the volume of the mixture.

        Parameters
        ----------
        moles : ndarray or list [float]
            Moles of each substance
        temperature : float
            System temperature [K]
        pressure : float
            System pressure [Pa]

        Returns
        -------
        float
            Volume of the mixture [m³]
        """
        mol_fractions = self.mol_fractions(moles)

        pure_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )
        return np.dot(pure_volumes, mol_fractions)

    def mix_heat_capacity(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Calculate heat capacity of th mixture.

        Parameters
        ----------
        moles : ndarray or list [float]
            Moles of each substance
        temperature : float
            System temperature [K]

        Returns
        -------
        float
            Heat capacity of the mixture [J/K]
        """
        mol_fractions = self.mol_fractions(moles)

        pure_cp = np.array(
            [
                substance.heat_capacity_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )
        mix_cp = np.dot(mol_fractions, pure_cp)
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

    def formation_enthalpies_correction(
        self, temperature: float, pressure: float
    ):
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation
        enthalpies of the pure substances from 298.15 K and 101325 Pa to
        the given temperature and pressure using Kirchhoff's equation.
        If the substance is a solid at a temperature > 298.15 K, this
        method also takes it into account.

        Parameters
        ----------
        temperature : float
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : float
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        correction_enthalpies : ndarray [float]
            Formation enthalpies correction for each substance (J/mol/K)
        """
        correction_enthalpies = np.array([])
        for substance in self.substances:
            if substance.normal_melting_point > 298.15:
                dhs = substance.heat_capacity_solid_dt_integral(
                    298.15, substance.normal_melting_point, pressure
                )
                dhf = substance.fusion_enthalpy(substance.normal_melting_point)
                dhl = substance.heat_capacity_liquid_dt_integral(
                    substance.normal_melting_point, temperature, pressure
                )

                correction_enthalpies = np.append(
                    correction_enthalpies, dhs + dhf + dhl
                )
            else:
                correction_enthalpies = np.append(
                    correction_enthalpies,
                    substance.heat_capacity_liquid_dt_integral(
                        298.15, temperature, pressure
                    ),
                )
        return correction_enthalpies

    def mixture_viscosity(
        self, temperature: float, pressure: float, moles: list
    ):
        """
        Evaluate viscosity of the mixture.

        Parameters
        ----------
        temperature : float
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : float
            Pressure at which formation enthalpies are to be calculated. [Pa]
        moles: list
        |   List of moles substance in the mixture
        Returns
        -------
         mixture_viscosity: float
            Viscosity of the mixture
        """
        raise NotImplementedError()
