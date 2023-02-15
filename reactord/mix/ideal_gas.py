"""Ideal gas Module."""
from typing import List

import chemicals

import numpy as np

from reactord.substance import Substance

from .abstract_mix import AbstractMix


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

    def __init__(self, substance_list: List[Substance]) -> None:
        self.substances = substance_list

    def volume(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> float:
        r"""Return the molar volume of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        .. math::
            v = \frac {R T} {P}

        | :math:`v`: mix's molar volume.
        | :math:`R`: Ideal gas constant = 8.31446261815324
          :math:`\frac {m^3 Pa} {K mol}`
        | :math:`T`: temperature.
        | :math:`P`: pressure.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        float
            Mixture's molar volume. [m³/mol]
        """
        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = r * np.divide(temperature, pressure)
        return volume

    def mix_heat_capacity(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> float:
        r"""Calculate the mixture's heat capacity [J/mol].

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        .. math::
            C_{p_{mix}} = \sum_{i=0}^{N} z_i C_{p_i}

        | :math:`C_{p_{mix}}`: mix's heat capacity.
        | :math:`N`: total number of substances in the mixture.
        | :math:`C_{p_i}`: ideal gas heat capacity of the mix's
          :math:`i`-th substance.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        float
            Mixture's heat capacity. [J/mol]


        Requires
        --------
            heat_capacity_gas defined on each mix's Substance.
        """
        pure_cp = np.array(
            [
                substance.heat_capacity_gas(temperature, pressure)
                for substance in self.substances
            ]
        )

        # The next line is equal to a column wise dot product of the two arrays
        mix_cp = np.multiply(mole_fractions, pure_cp).sum(axis=0)
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

    def formation_enthalpies(self, temperature: float, pressure: float):
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation enthalpies
        of the pure substances from 298.15 K and 101325 Pa to the given
        temperature and pressure.

        .. math::
           \Delta H_{f T_{(g)}} = \Delta H_{f 298.15_{(g)}}^0

        | :math:`C_{p_{mix}}`: mix's heat capacity.
        | :math:`N`: total number of substances in the mixture.
        | :math:`C_{p_i}`: ideal gas heat capacity of the mix's
          :math:`i`-th substance.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.

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
