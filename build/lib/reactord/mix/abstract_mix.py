"""Abstract_mix Module."""
from abc import ABCMeta, abstractmethod
from typing import Callable, List, Union

import numpy as np
from numpy.typing import NDArray

from reactord.mix.viscosity_mixing_rules.grunbergnissan import grunberg_nissan
from reactord.mix.viscosity_mixing_rules.herningzipperer import (
    herning_zipperer,
)
from reactord.mix.viscosity_mixing_rules.linearmix import linear
from reactord.substance.substance import Substance


class AbstractMix(metaclass=ABCMeta):
    """Mixture objects abstract class.

    Parameters
    ----------
    metaclass : AbstractMix
        Mixture objects interface.
    phase_nature : str
        Liquid or gas.
    substance_list : List[Substance]
        List of substances.
    viscosity_mixing_rule : str
        Viscosity mixing rule method. Options available: "linear",
        "grunberg_nissan", "herning_zipperer".

    Attributes
    ----------
    names : NDArray[str]
       Names of the mixture's substances.
    molecular_weights : NDArray[np.float64]
        The molecular weights of the mixture's substances. [g/mol]
    normal_boiling_points : NDArray[np.float64]
        The normal boiling points of the mixture's substances. [K]
    normal_melting_points : NDArray[np.float64]
        The normal melting points of the mixture's substances. [K]
    critical_temperatures : NDArray[np.float64]
        The critical temperatures of the mixture's substances. [K].
    critical_pressures : NDArray[np.float64]
        The critical pressures of the mixture's substances. [Pa]
    acentric_factors : NDArray[np.float64]
        The acentric factors of the mixture's substances.
    formation_enthalpies : NDArray[np.float64]
        Standard state molar enthalpies of formation of the mixture's
        substances. [J/mol]
    formation_enthalpies_ig : NDArray[np.float64]
        Ideal-gas molar enthalpies of formation of the mixture's substances.
        [J/mol]
    """

    substances: List[Substance] = []
    phase_nature: str = "gas"
    viscosity_mixing_rule: str = "linear"
    _viscosity_mixing_rule_function: Callable = linear

    def __init__(
        self,
        substance_list: List[Substance],
        phase_nature: str,
        viscosity_mixing_rule: str,
    ) -> None:
        self.substances = substance_list
        self.phase_nature = phase_nature.lower()
        self.viscosity_mixing_rule = viscosity_mixing_rule.lower()

        # Check if phase nature is correct
        if self.phase_nature not in ["gas", "liquid"]:
            raise ValueError(
                f"{self.phase_nature} is not a valid phase nature"
            )

        # Choose the function for viscosity mixing rule.
        # TODO: TAL VEZ ESTO DEBERIA SER ALGO INYECTABLE SINO ESTAS LIMITADO AL
        # HARDCODEO DEL SOURCE DE AbstractMix
        if self.viscosity_mixing_rule == "grunberg_nissan":
            self._viscosity_mixing_rule_function = grunberg_nissan
        elif self.viscosity_mixing_rule == "herning_zipperer":
            self._viscosity_mixing_rule_function = herning_zipperer
        elif self.viscosity_mixing_rule == "linear":
            self._viscosity_mixing_rule_function = linear
        else:
            raise ValueError(
                f"""{self.viscosity_mixing_rule} is not a valid viscosity
                mixing rule."""
            )

        # =====================================================================
        # Put substances attributes in arrays.
        # =====================================================================
        self.names = np.array([])
        self.molecular_weights = np.array([])
        self.normal_boiling_points = np.array([])
        self.normal_melting_points = np.array([])
        self.critical_temperatures = np.array([])
        self.critical_pressures = np.array([])
        self.acentric_factors = np.array([])
        self.formation_enthalpies = np.array([])
        self.formation_enthalpies_ig = np.array([])

        for substance in self.substances:
            self.names = np.append(self.names, substance.name)
            self.molecular_weights = np.append(
                self.molecular_weights, substance.molecular_weight
            )
            self.normal_boiling_points = np.append(
                self.normal_boiling_points, substance.normal_boiling_point
            )
            self.normal_melting_points = np.append(
                self.normal_melting_points, substance.normal_melting_point
            )
            self.critical_temperatures = np.append(
                self.critical_temperatures, substance.critical_temperature
            )
            self.critical_pressures = np.append(
                self.critical_pressures, substance.critical_pressure
            )
            self.acentric_factors = np.append(
                self.acentric_factors, substance.acentric_factor
            )
            self.formation_enthalpies = np.append(
                self.formation_enthalpies, substance.formation_enthalpy
            )
            self.formation_enthalpies_ig = np.append(
                self.formation_enthalpies_ig, substance.formation_enthalpy_ig
            )

        # Check for duplicates:
        if np.size(np.unique(self.names)) != np.size(self.names):
            raise ValueError("All mix's substances must have different names.")

    # =========================================================================
    # Mixtures' common methods
    # =========================================================================
    def mole_fractions(
        self, moles: Union[float, NDArray[np.float64]]
    ) -> Union[float, NDArray[np.float64]]:
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
        moles: Union[float, NDArray[np.float64]]
            Mole fractions of each substance specified in the same order as the
            mix's substances order.

        Returns
        -------
        Union[float, NDArray[np.float64]]
            The array of the molar fractions of the mix's substances. Each
            column of the matrix represents each mixture and each row represent
            each substance's mole fraction.
        """
        total_moles = np.sum(moles, axis=0)
        mole_fractions = np.divide(moles, total_moles)

        return mole_fractions

    def mix_molecular_weight(
        self, mole_fractions: Union[float, NDArray[np.float64]]
    ) -> Union[float, NDArray[np.float64]]:
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
        mole_fractions : Union[float, NDArray[np.float64]]
            Mole fractions of each substance specified in the same order as the
            mix's substances order.

        Returns
        -------
        Union[float, NDArray[np.float64]]
            The molecular weight of the mixture. [kg/kmol]


        Requires
        --------
            molecular_weight defined on each mix's Substance.
        """
        return np.dot(self.molecular_weights, mole_fractions)

    def molar_density(
        self,
        mole_fractions: Union[float, NDArray[np.float64]],
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
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
        mole_fractions : Union[float, NDArray[np.float64]]
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's molar density. [mol/m3]
        """
        molar_volume = self.volume(mole_fractions, temperature, pressure)

        return np.divide(1, molar_volume)

    def mass_density(
        self,
        mole_fractions: Union[float, NDArray[np.float64]],
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
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
        mole_fractions : Union[float, NDArray[np.float64]]
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's mass density. [kg/m3]


        Requires
        --------
            molecular_weight defined on each mix's Substance.
        """
        mass_density = (
            self.molar_density(mole_fractions, temperature, pressure)
            * self.mix_molecular_weight(mole_fractions)
            / 1000
        )

        return mass_density

    def mix_viscosity(
        self,
        mole_fractions: NDArray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
        """Evaluate the viscosity of the mixture with the chosen mixing rule.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : NDArray
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's viscosity. [Pa s]
        """
        mixture_viscosities = self._viscosity_mixing_rule_function(
            self, mole_fractions, temperature, pressure
        )

        return mixture_viscosities

    # Compositional measure functions
    def partial_pressures(
        self,
        mole_fractions: NDArray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
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
        mole_fractions : NDArray
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            An array that contains the partial pressures of mix's substances.
            Each column of the matrix represents each mixture and each row
            represent each substance's partial pressure. [Pa]
        """
        partial_pressures = np.multiply(mole_fractions, pressure)

        return partial_pressures

    def concentrations(
        self,
        mole_fractions: NDArray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
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
        mole_fractions : NDArray
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            An array that contains the concentrations of the mixture's
            substances. [mol/m³]
        """
        molar_densities = self.molar_density(
            mole_fractions, temperature, pressure
        )

        return np.multiply(molar_densities, mole_fractions)

    def __len__(self) -> int:
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
        self,
        mole_fractions: NDArray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
        """Return the molar volume of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : NDArray
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's molar volume. [m³/mol]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def mix_heat_capacity(
        self,
        mole_fractions: NDArray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
        """Calculate the mixture's heat capacity.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        Parameters
        ----------
        mole_fractions : NDArray
            Mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's heat capacity. [J/mol]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def formation_enthalpies_correction(
        self,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ) -> Union[float, NDArray[np.float64]]:
        """Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation enthalpies
        of the pure substances from 298.15 K and 101325 Pa to the given
        temperature and pressure. Multiple temperatures and pressures may be
        specified following de NumPy broadcasting rules.

        Parameters
        ----------
        temperature : Union[float, NDArray[np.float64]]
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : Union[float, NDArray[np.float64]]
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Formation enthalpies correction for each mix's substance. [J/mol/K]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @abstractmethod
    def get_formation_enthalpies(self) -> Union[float, NDArray[np.float64]]:
        """Calculate the formation enthalpy of the mixture's substances.

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Formation enthalpies of the mixture's substances. [J/mol]

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")
