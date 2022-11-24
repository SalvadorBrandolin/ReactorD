from abc import ABCMeta, abstractmethod
from typing import List

import numpy as np

from reactord.substance import Substance


class AbstractMix(metaclass=ABCMeta):
    """Mixture object abstract class.

    Parameters
    ----------
    substance_list : ndarray or List[Substance objects]
        list or array of Substance objects."""

    substances: List[Substance] = []

    # ==================================================================
    # Mixtures common methods
    # ==================================================================
    def mol_fracations(self, moles: List[float]):
        """method that calculates the molar fractions of the mixture

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance

        Returns
        -------
        ndarray
            array that contains the molar fractions of mixture's
            substances
        """
        total_moles = np.sum(moles, axis=0)
        zi = np.divide(moles, total_moles)

        return zi

    def partial_pressures(
        self, moles: List[float], temperature: float, pressure: float
    ):
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

    def __len__(self):
        return len(self.substances)

    def __str__(self):
        string = (
            f"The mixture contains the following"
            f" {len(self.substances)} components:\n"
        )
        for i, substance in enumerate(self.substances):
            string = string + substance.name.capitalize() + "\n"
        return string

    # ==================================================================
    # Mixtures specifics methods
    # ==================================================================

    @abstractmethod
    def concentrations(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Concentrations of the mixtures substances at the given moles
        of each compound, temperature and pressure.

        Parameters:
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]

        Returns
        -------
        ndarray [float]
            ndarray that contains the concentrations of the mixture's
            substances [mol/m^3]
        """
        raise NotImplementedError()

    @abstractmethod
    def volume(self):
        """Method that returns the volume of the mixture.

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
        float
            volume of the mixture [m^3]
        """
        raise NotImplementedError()

    @abstractmethod
    def mix_heat_capacity(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Method that returns the heat capacity of the mixture.

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]
        """
        raise NotImplementedError()

    @abstractmethod
    def _formation_enthalpies_set(self):
        """Method that obtains or calculates the formation enthalpy of the
        mixture's substances.
        """
        raise NotImplementedError()

    @abstractmethod
    def formation_enthalpies_correction(
        self,
        temperature: float,
        pressure: float,
    ):
        """Method that calculates the correction therm for the formation
        enthalpy of the pure substances from 298.15 K and 100000 Pa to
        temperature and pressure

        Parameters
        ----------
        temperature : float
            Correction temperature for the formation enthalpies. [K]
        pressure : float
            Correction pressure for the formation enthalpies. [Pa]

        """
        raise NotImplementedError()
