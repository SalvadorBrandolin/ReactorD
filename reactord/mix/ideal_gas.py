"""Ideal gas Module."""
from typing import List

import chemicals

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
        moles : ndarray or list [float]
            Moles of each substance
        temperature : float
            System temperature [K]
        pressure : float
            System pressure [Pa]

        Returns
        -------
        ndarray or list [float]
            Concentration of each substance [mol/m³]
        """
        mol_fractions = self.mol_fractions(moles)

        r = 8.31446261815324  # m³⋅Pa/K/mol
        density = pressure / (r * temperature)
        return np.multiply(mol_fractions, density)

    def volume(
        self, moles: List[float], temperature: float, pressure: float
    ) -> float:
        """Calculate the volume of the mixture.

        Parameters
        ----------
        moles : ndarray or list[float]
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
        total_moles = np.sum(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = total_moles * r * temperature / pressure
        return volume

    def mix_heat_capacity(
        self, moles: List[float], temperature: float, pressure: float
    ) -> float:
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
                substance.heat_capacity_gas(temperature, pressure)
                for substance in self.substances
            ]
        )
        mix_cp = np.dot(mol_fractions, pure_cp)
        return mix_cp

    def _formation_enthalpies_set(self):
        """Return the ideal gas formation enthalpies in a ordered ndarray.

        Method that read the ideal gas formation enthalpies of the mix
        class and returns them in a ordered ndarray.

        Returns
        -------
        ndarray [float]
            Ideal gas formation enthalpies of each substance [J/mol/K]
        """
        enthalpies = np.array([])

        for substance in self.substances:
            enthalpies = np.append(enthalpies, substance.formation_enthalpy_ig)

        return enthalpies

    def formation_enthalpies_correction(
        self, temperature: float, pressure: float
    ):
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation
        enthalpies of the pure substances from 298.15 K and 101325 Pa to
        the given temperature and pressure using Kirchhoff's equation.

        Parameters
        ----------
        temperature : float
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : float
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        correction_enthalpies : ndarray [float]
            Formation enthalpies of each substance (J/mol/K)
        """
        correction_enthalpies = np.array([])

        for substance in self.substances:
            correction_enthalpies = np.append(
                correction_enthalpies,
                substance.heat_capacity_gas_dt_integral(
                    298.15, temperature, pressure
                ),
            )

        return correction_enthalpies

    def mixture_viscosity(
        self,
        moles: List[float],
        temperature: float,
        pressure: float,
    ) -> float:
        """
        Evaluate the viscosity of the mixture.

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
        mol_fractions = self.mol_fractions(moles)
        viscosity_pure = []
        molecular_weight = []
        for substance in self.substances:
            viscosity_pure.append(
                substance.viscosity_gas(temperature, pressure)
            )
            molecular_weight.append(substance.molecular_weight)

        return chemicals.viscosity.Herning_Zipperer(
            mol_fractions, viscosity_pure, molecular_weight
        )
