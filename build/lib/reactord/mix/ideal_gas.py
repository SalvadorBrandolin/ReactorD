"""Ideal gas Module."""
from typing import List, Union

import numpy as np
from numpy.typing import NDArray

from reactord.substance import Substance

from scipy.constants import R

from .abstract_mix import AbstractMix


class IdealGas(AbstractMix):
    """IdealGas object class.

    This class should be used when the mixture is considered an ideal gas.

    Parameters
    ----------
    substance_list : list[Substance]
        list of substance objects
    viscosity_mixing_rule : str, optional
        Viscosity mixing rule method. Options available: "linear",
        "grunberg_nissan", "herning_zipperer"., by default "herning_zipperer"

    Attributes
    ----------
    phase_nature : str
        gas.
    substances : List [Substance]
        List of substances.
    viscosity_mixing_rule : str
        Viscosity mixing rule method.
    """

    def __init__(
        self,
        substance_list: List[Substance],
        viscosity_mixing_rule: str = "herning_zipperer",
    ) -> None:
        super().__init__(
            substance_list=substance_list,
            phase_nature="gas",
            viscosity_mixing_rule=viscosity_mixing_rule,
        )

    def volume(
        self,
        mole_fractions: np.ndarray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
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
        mole_fractions : NDArray[np.float64]
            mole fractions of each substance specified in the same order as the
            mix's substances order.
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Mixture's molar volume. [m³/mol]
        """
        volume = R * np.divide(temperature, pressure)
        return volume

    def mix_heat_capacity(
        self,
        mole_fractions: np.ndarray,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
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
        temperature: Union[float, NDArray[np.float64]]
            Temperature. [K]
        pressure: Union[float, NDArray[np.float64]]
            Pressure. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
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

        # The next line is equal to a column-wise dot product of the two arrays
        mix_cp = np.multiply(mole_fractions, pure_cp).sum(axis=0)
        return mix_cp

    def formation_enthalpies_correction(
        self,
        temperature: Union[float, NDArray[np.float64]],
        pressure: Union[float, NDArray[np.float64]],
    ):
        r"""Calculate the correction term for the formation enthalpy.

        Method that calculates the correction term for the formation enthalpies
        of the pure substances from 298.15 K and 101325 Pa to the given
        temperature and pressure.

        .. math::
           \Delta H_{f T_{(g)}} = \Delta H_{f 298.15_{(g)}}^0 + \sum_{i=0}^{N}
           \int_{298.15}^{T} v_i C_{p_i} dT

           CT = \sum_{i=0}^{N} \int_{298.15}^{T} v_i C_{p_i} dT

        | :math:`\Delta H_{f T_{(g)}}`: reaction enthalpy.
        | :math:`\Delta H_{f 298.15_{(g)}}^0`: ideal gas standard reaction
          enthalpy.
        | :math:`N`: total number of substances in the mixture.
        | :math:`C_{p_i}`: ideal gas heat capacity of the mix's
          :math:`i`-th substance.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
        | :math:`v_i`: stoichiometric coefficient of the mix's :math:`i`-th
          substance.
        | :math:`CT`: correction term.

        Parameters
        ----------
        temperature : Union[float, NDArray[np.float64]]
            Temperature at which formation enthalpies are to be calculated. [K]
        pressure : Union[float, NDArray[np.float64]]
            Pressure at which formation enthalpies are to be calculated. [Pa]

        Returns
        -------
        Union[float, NDArray[np.float64]]
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

    def get_formation_enthalpies(self):
        """Return the ideal gas formation enthalpies in an ordered ndarray.

        Method that read the ideal gas formation enthalpies of the mix
        class and returns them in an ordered ndarray.

        Returns
        -------
        Union[float, NDArray[np.float64]]
            Ideal gas formation enthalpies of each substance [J/mol/K]
        """
        return self.formation_enthalpies_ig
