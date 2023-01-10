"""Abstract_mix Module."""
from abc import ABCMeta, abstractmethod
from typing import List

import numpy as np

from reactord.substance import Substance


class AbstractMix(metaclass=ABCMeta):
    """Mixture object abstract class.

    Parameters
    ----------
    metaclass : AbstractMix, optional
        Mixture objects interface

    Raises
    ------
    NotImplementedError
        concentrations abstract method not implemented
    NotImplementedError
        volume abstract method not implemented
    NotImplementedError
        mix_heat_capacity abstract method not implemented
    NotImplementedError
        _formation_enthalpies_set abstract method not implemented
    NotImplementedError
        formation_enthalpies_correction abstract method not implemented
    """

    substances: List[Substance] = []

    # ==================================================================
    # Mixtures common methods
    # ==================================================================

    def mol_fractions(self, moles: List[float]):
        """Calculate the molar fractions of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        row represent each substance and each colum represent each mixture
        composition.

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance specified in the same order as the
            mix substances order.

        Returns
        -------
        mol_fractions : ndarray of shape (moles,)
        Array of the molar fractions of mixture's substances
        """
        total_moles = np.sum(moles, axis=0)
        mol_fractions = np.divide(moles, total_moles)

        return mol_fractions

    def partial_pressures(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Calculate the partial pressures of the mixture's components.

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
        partial_pressures: ndarray of shape(moles,) [Pa]
            array that contains the partial pressures of mixture's
            substances
        """
        mol_fractions = self.mol_fractions(moles)
        partial_pressures = np.multiply(mol_fractions, pressure)

        return partial_pressures

    def __len__(self):
        """Return the number of substances in the mixture.

        Redifine the magic method __len__
        Returns
        -------
        int
            Number of substances in the mixture
        """
        return len(self.substances)

    def __str__(self):
        """Return a text with information about the components of the mixture.

        Returns
        -------
        str
            Text with information about the components of the mixture.
        """
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
        """Concentrations of the mixture's substances.

        Concentrations of the mixture's substances given moles
        of each compound, temperature and pressure.

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
        ndarray [float] of shape (moles,)
            ndarray that contains the concentrations of the mixture's
            substances [mol/m³]
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def volume(self):
        """Return the volume of the mixture.

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
            volume of the mixture [m³]
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def mix_heat_capacity(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Return the heat capacity of the mixture.

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
            heat capacity of the mixture [J/K]

        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def formation_enthalpies_correction(
        self,
        temperature: float,
        pressure: float,
    ):
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation
        enthalpies of the pure substances from 298.15 K and 101325 Pa to
        the given temperature and pressure

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
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def _formation_enthalpies_set(self):
        """Calculate the formation enthalpy of the mixture's substances."""
        raise NotImplementedError("Abstract method not implemented")
