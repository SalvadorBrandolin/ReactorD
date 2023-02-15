"""Abstract_mix Module."""
from abc import ABCMeta, abstractmethod
from typing import List

import numpy as np

from reactord.substance import Substance


class AbstractMix(metaclass=ABCMeta):
    """Mixture objects abstract class.

    Parameters
    ----------
    metaclass : AbstractMix
        Mixture objects interface.

    Attributes
    ----------
    substances : list [Substance]
        List of substances.
    viscosity_mixing_rule : str
        Viscosity mixing rule method. Method available:

    Raises
    ------
    NotImplementedError
        concentrations abstract method not implemented.
    NotImplementedError
        volume abstract method not implemented.
    NotImplementedError
        mix_heat_capacity abstract method not implemented.
    NotImplementedError
        formation_enthalpies_correction abstract method not implemented.
    NotImplementedError
        _formation_enthalpies_set abstract method not implemented.
    NotImplementedError
        mix_viscosity abstract method not implemented.
    """

    substances: List[Substance] = []
    viscosity_mixing_rule : str = "linear"

    # =========================================================================
    # Mixtures' common methods
    # =========================================================================
    def mole_fractions(self, moles: np.ndarray) -> np.ndarray:
        r"""Calculate the mole fractions of the mixture's substances.

        Multiple mixture compositions can be specified by a moles matrix. Each 
        column of the matrix represents each mixture and each row represent 
        each substance's mole fraction.

        .. math::
            z_i = \frac {n_i} {\sum_{i=0}^{N} n_i}

        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
        | :math:`n_i`: moles of the mix's :math:`i`-th substance.
        | :math:`N`: total number of substances in the mixture.

        Parameters
        ----------
        moles: np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.

        Returns
        -------
        ndarray
            The array of the molar fractions of the mix's substances. Each 
            column of the matrix represents each mixture and each row represent
            each substance's mole fraction.
        """
        total_moles = np.sum(moles, axis=0)
        mole_fractions = np.divide(moles, total_moles)

        return mole_fractions

    def mix_molecular_weight(self, mole_fractions: np.ndarray) -> float:
        r"""Calculate the molecular weight of the mixture.

        Multiple mixture compositions can be specified by a mole_fractions
        matrix. Each column of the matrix represents each mixture and each row 
        represent each substance's mole fraction. In this case, the return
        is a 1D array with the molecular weight of each mix.

        .. math::
            M_{mix} = \sum_{i=0}^{N} M_i z_i

        | :math:`M_{mix}`: molecular weight of the mixture.
        | :math:`N`: total number of substances in the mixture.
        | :math:`M_i`: molecular weight of the mix's :math:`i`-th substance.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the
            mix's substances order.

        Returns
        -------
        float
            The molecular weight of the mixture. [kg/kmol]


        Requires
        --------
            molecular_weight defined on each mix's Substance.
        """
        pure_molecular_weights = np.array(
            [substance.molecular_weight for substance in self.substances]
        )

        return np.dot(pure_molecular_weights, mole_fractions)

    def molar_density(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> float:
        r"""Calculate the mixture's molar density.

        Multiple mixture compositions can be specified by a mole_fractions
        matrix. Each column of the matrix represents each mixture and each row 
        represent each substance's mole fraction. Also with temperature and 
        pressure vectors following de NumPy broadcasting rules.

        .. math::
            \rho_{mix} = \frac {1} {v_{(\vec{z}, T, P)}}

        | :math:`\rho_{mix}`: mix's molar density.
        | :math:`v`: mix's molar volume.
        | :math:`T`: temperature.
        | :math:`P`: pressure.
        | :math:`\vec{z}`: substances' mole fractions.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        float
            Mixture's molar density. [mol/m3]
        """
        molar_volume = self.volume(mole_fractions, temperature, pressure)

        return np.divide(1, molar_volume)

    def mass_density(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> float:
        r"""Calculate the mixture's mass density.

        Multiple mixture compositions can be specified by a mole_fractions
        matrix. Each column of the matrix represents each mixture and each row 
        represent each substance's mole fraction. Also with temperature and 
        pressure vectors following de NumPy broadcasting rules.

        .. math::
            \rho_{mix}^{mass} = \rho_{mix}^{mol} M_{mix}

        | :math:`\rho_{mix}^{mass}`: mix's masic density.
        | :math:`\rho_{mix}^{mol}`: mix's molar density.
        | :math:`M_{mix}`: molecular weight of the mixture.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        float
            Mixture's mass density. [kg/m3]


        Requires
        --------
            molecular_weight defined on each mix's Substance.
        """
        mass_density = self.molar_density(
            mole_fractions, temperature, pressure
        ) * self.mix_molecular_weight(mole_fractions)

        return mass_density

    # Compositional measure functions
    def partial_pressures(
        self,
        mole_fractions: np.ndarray,
        temperature: float,
        pressure: float,
    ) -> np.ndarray:
        """Calculate the partial pressures of the mixture's substances.

        Multiple mixture compositions can be specified by a mole_fractions
        matrix. Each column of the matrix represents each mixture and each row 
        represent each substance's mole fraction. Also with temperature and 
        pressure vectors following de NumPy broadcasting rules.

        .. math::
            P_i = P z_i

        | :math:`P_i`: partial pressure of the :math:`i`-th mix's substance.
        | :math:`P`: total system's pressure.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        ndarray [float]
            An array that contains the partial pressures of mix's substances.
            Each column of the matrix represents each mixture and each row 
            represent each substance's partial pressure. [Pa]
        """
        partial_pressures = np.multiply(mole_fractions, pressure)

        return partial_pressures

    def concentrations(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> np.ndarray:
        r"""Concentrations of the mixture's substances.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure 
        vectors following de NumPy broadcasting rules.

        .. math::
            C_i = z_i \rho_{mix}

        | :math:`C_i`: concentration of the mix's :math:`i`-th substance.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
        | :math:`\rho_{mix}`: mix's molar density.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Returns
        -------
        ndarray [float]
            An array that contains the concentrations of the mixture's
            substances. [mol/m³]
        """
        molar_volumes = self.volume(mole_fractions, temperature, pressure)

        return np.multiply(mole_fractions, molar_volumes)

    def __len__(self):
        """Return the number of substances in the mixture.

        Redifine the magic method __len__

        Returns
        -------
        int
            The number of substances in the mixture.
        """
        return len(self.substances)

    # =========================================================================
    # Mixtures specifics methods
    # =========================================================================
    @abstractmethod
    def volume(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> None:
        """Return the molar volume of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure 
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def mix_heat_capacity(
        self, mole_fractions: np.ndarray, temperature: float, pressure: float
    ) -> None:
        """Calculate the mixture's heat capacity.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure 
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

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
        temperature and pressure. Multiple temperatures and pressures may be 
        specified following de NumPy broadcasting rules.

        Parameters
        ----------
        temperature : float
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : float
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        ndarray [float]
            Formation enthalpies correction for each mix's substance. [J/mol/K]

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
    def mix_viscosity(
        self,
        mole_fractions: np.ndarray,
        temperature: float,
        pressure: float,
    ) -> None:
        """
        Evaluate the viscosity of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure 
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : np.ndarray [float]
            Mole fractions of each substance specified in the same order as the 
            mix's substances order.
        temperature: float
            Temperature. [K]
        pressure: float
            Pressure. [Pa]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")
