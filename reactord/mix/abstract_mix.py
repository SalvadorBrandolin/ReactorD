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
        Mixture objects interface.

    Raises
    ------
    NotImplementedError
        concentrations abstract method not implemented.
    NotImplementedError
        volume abstract method not implemented.
    NotImplementedError
        mix_heat_capacity abstract method not implemented.
    NotImplementedError
        _formation_enthalpies_set abstract method not implemented.
    NotImplementedError
        formation_enthalpies_correction abstract method not implemented.
    """

    substances: List[Substance] = []

    # =========================================================================
    # Mixtures common methods
    # =========================================================================

    def mol_fractions(self, moles: List[float]) -> np.ndarray:
        """Calculate the molar fractions of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        row represent each substance and each colum represent each mixture
        composition.

        Parameters
        ----------
        moles: ndarray or list [float]
            Moles of each substance specified in the same order as the mix
            substances order.

        Returns
        -------
        mol_fractions : ndarray [float]
            Array of the molar fractions of mixture's substances.
        """
        total_moles = np.sum(moles, axis=0)
        mol_fractions = np.divide(moles, total_moles)

        return mol_fractions

    def mixture_molecular_weight(self, moles: List[float]):
        """Calculate the molecular weight of the mixture.

        Parameters
        ----------
        moles : List[float]
            moles of each substance specified in the same order as the
            mix substances order.

        Returns
        -------
        mixture_molecular_weight: ndarray [g/mol]
            molecular weight of the mixture calculated based on molar fractions
        """
        pure_molecular_weights = [
            substance.molecular_weight for substance in self.substances
        ]

        molar_fractions = self.mol_fractions(moles)

        return np.dot(pure_molecular_weights, molar_fractions)

    def partial_pressures(
        self, moles: List[float], temperature: float, pressure: float
    ) -> np.ndarray:
        """Calculate the partial pressures of the mixture's components.

        Parameters
        ----------
        moles: ndarray or list [float]
            Moles of each substance.
        temperature: float
            Temperature. [K]
        pressure: float
            Total Pressure. [Pa]

        Returns
        -------
        partial_pressures: ndarray [float]
            Array that contains the partial pressures of mixture's substances.
            [Pa]
        """
        mol_fractions = self.mol_fractions(moles)
        partial_pressures = np.multiply(mol_fractions, pressure)

        return partial_pressures

    def molar_density(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Calculate mixture molar density [mol/m3].

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
        molar_density: float
        """
        total_moles = np.sum(moles)
        return total_moles / self.volume(moles, temperature, pressure)

    def mass_density(
        self, moles: List[float], temperature: float, pressure: float
    ):
        """Return density in [kg/m3].

        Calculate mixture mass density.

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
        mass_density: float
        """
        mass_density = (
            self.molar_density(moles, temperature, pressure)
            * self.mixture_molecular_weight(moles)
            / 1000
        )

        return mass_density

    def __len__(self):
        """Return the number of substances in the mixture.

        Redifine the magic method __len__

        Returns
        -------
        int
            Number of substances in the mixture.
        """
        return len(self.substances)

    def __str__(self) -> str:
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

    # =========================================================================
    # Mixtures specifics methods
    # =========================================================================
    @abstractmethod
    def concentrations(
        self, moles: List[float], temperature: float, pressure: float
    ) -> None:
        """Concentrations of the mixture's substances.

        Concentrations of the mixture's substances given moles of each
        compound, temperature and pressure.

        Parameters
        ----------
        moles: ndarray or list [float]
            Moles of each substance.
        temperature: float
            Temperature. [K]
        pressure: float
            Total Pressure. [Pa]

        Returns
        -------
        ndarray [float]
            ndarray that contains the concentrations of the mixture's
            substances. [mol/m³]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def volume(
        self, moles: List[float], temperature: float, pressure: float
    ) -> None:
        """Return the volume of the mixture.

        Parameters
        ----------
        moles: ndarray or list [float]
            Moles of each substance.
        temperature: float
            Temperature. [K]
        pressure: float
            Total Pressure. [Pa]

        Returns
        -------
        float
            volume of the mixture. [m³]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def mix_heat_capacity(
        self, moles: List[float], temperature: float, pressure: float
    ) -> None:
        """Return the heat capacity of the mixture.

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance.
        temperature: float
            Temperature. [K]
        pressure: float
            Total Pressure. [Pa]

        Returns
        -------
        float
            heat capacity of the mixture. [J/K]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def formation_enthalpies_correction(
        self,
        temperature: float,
        pressure: float,
    ) -> None:
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation enthalpies
        of the pure substances from 298.15 K and 101325 Pa to the given
        temperature and pressure.

        Parameters
        ----------
        temperature : float
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : float
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        correction_enthalpies : ndarray [float]
            Formation enthalpies correction for each substance. [J/mol/K]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def _formation_enthalpies_set(self) -> None:
        """Calculate the formation enthalpy of the mixture's substances."""
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def mixture_viscosity(
        self,
        moles: list,
        temperature: float,
        pressure: float,
    ):
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
        raise NotImplementedError()
