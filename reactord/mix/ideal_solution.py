"""Ideal solution Module."""
from typing import List

import numpy as np

from reactord.substance import Substance

from .abstract_mix import AbstractMix


class IdealSolution(AbstractMix):
    """IdealSolution object class.

    This class should be used when the mixture is considered an ideal solution.

    Parameters
    ----------
    substance_list : list [float]
        list of substance objects
    viscosity_mixing_rule : str, optional
        Viscosity mixing rule method. Options available: "linear",
        "grunberg_nissan", "herning_zipperer"., by default "grunberg_nissan"

    Attributes
    ----------
    phase_nature : str
        liquid.
    substances : List [Substance]
        List of substances.
    viscosity_mixing_rule : str
        Viscosity mixing rule method.
    """

    def __init__(
        self,
        substance_list: List[Substance],
        viscosity_mixing_rule: str = "grunberg_nissan",
    ) -> None:
        super().__init__(
            substance_list=substance_list,
            phase_nature="liquid",
            viscosity_mixing_rule=viscosity_mixing_rule,
        )

    def volume(
        self, mole_fractions: List[float], temperature: float, pressure: float
    ) -> float:
        r"""Return the molar volume of the mixture.

        Multiple mixture compositions can be specified by a moles matrix. Each
        column of the matrix represents each mixture and each row represents
        each substance's mole fractions. Also with temperature and pressure
        vectors following de NumPy broadcasting rules.

        .. math::
            v = \sum z_i {v_i}_{(T,P)}

        | :math:`v`: mix's molar volume.
        | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
        | :math:`{v_i}_{(T,P)`: liquid molar volume of the mix's :math:`i`-th
          substance.
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


        Requires
        --------
            volume_liquid defined on each mix's Substance.
        """
        pure_volumes = np.array(
            [
                substance.volume_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )

        mix_volumes = np.multiply(pure_volumes, mole_fractions).sum(axis=0)

        return mix_volumes

    def mix_heat_capacity(
        self, mole_fractions: List[float], temperature: float, pressure: float
    ):
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
            heat_capacity_liquid defined on each mix's Substance.
        """
        pure_cp = np.array(
            [
                substance.heat_capacity_liquid(temperature, pressure)
                for substance in self.substances
            ]
        )

        # The next line is equal to a column-wise dot product of the two arrays
        mix_cp = np.multiply(mole_fractions, pure_cp).sum(axis=0)
        return mix_cp

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

    def get_formation_enthalpies(self):
        """Return the formation enthalpies in a ordered ndarray.

        Method that read the formation enthalpies of mix and returns
        them in a ordered ndarray.

        Returns
        -------
        ndarray [float]
            Formation enthalpies of each substance
        """
        return self.formation_enthalpies
