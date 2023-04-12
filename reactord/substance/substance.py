"""Substance module.

Class to define a substance for ReactorD library.
"""

from typing import Callable

from dill import dumps, loads

import numpy as np

from .symbolic import Symbolic
from .thermo_substance import thermo_substance_constructor


class Substance(Symbolic):
    """Substance object class.

    Class to define a substance object. Specific attributes definition will be
    required for the reactors, described in each reactor documentation. For
    example, an adiabatic reactor will require that substances define a heat
    capacity function, on the other hand, when using isothermic reactors this
    won't be necessary. Substance has the from_thermo_data_base alternative
    construction method. Substances objects can be saved as pickle files with
    the method to_pickle. substances can be loaded from a pickle file with the
    method from_pickle

    Example:

    >>> water = Substance.from_thermo_data_base('water')
    >>> water.to_pickle('my_water_file')
    >>> water_pickle = water.from_pickle('my_water_file)

    Parameters
    ----------
    name : str, optional
       Name of the substance, by default None.
    molecular_weight : float, optional
        The molecular weight of the substance [g/mol], by default None.
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None.
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None.
    critical_temperature : float, optional
        The critical temperature of the substance [K], by default None.
    critical_pressure : float, optional
        The critical pressure of the substance [Pa], by default None.
    acentric_factor : float, optional
        The acentric factor of the substance, by default None.
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None.
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None.
    vaporization_enthalpy : Callable, optional
        A function that receives a temperature and returns the vaporization
        enthalpy at temperature [J/mol], by default None.
    sublimation_enthalpy : Callable, optional
        A function that receives temperature and returns the sublimation
        enthalpy at temperature [J/mol], by default None.
    volume_solid : Callable, optional
        A function that receives temperature and pressure, and returns the
        molar volume of the solid at temperature [m³/mol], by default None.
    volume_liquid : Callable, optional
        A function that receives temperature and pressure, and returns the
        molar volume of liquid at temperature and pressure [m³/mol], by default
        None.
    heat_capacity_solid : Callable, optional
        A function that receives temperature and pressure, and returns the
        heat capacity of the solid at temperature and pressure [J/mol/K], by
        default None.
    heat_capacity_liquid : Callable, optional
        A function that receives temperature and pressure, and returns the
        heat capacity of the liquid at temperature and pressure [J/mol/K], by
        default None.
    heat_capacity_gas : Callable, optional
        A function that receives temperature and pressure, and returns the
        ideal gas heat capacity at temperature and pressure [J/mol/K], by
        default None.
    thermal_conductivity_liquid : Callable, optional
        A function that receives temperature and pressure, and returns the
        thermal conductivity of the liquid at temperature and pressure [W/m/K],
        by default None.
    thermal_conductivity_gas : Callable, optional
        A function that receives temperature and pressure, and returns the
        thermal conductivity of the gas at temperature and pressure [W/m/K], by
        default None.
    viscosity_liquid : Callable, optional
        A function that receives temperature and pressure, and returns the
        viscosity of the liquid at temperature and pressure [Pa*s], by default
        None.
    viscosity_gas : Callable, optional
        A function that receives temperature and pressure, and returns the
        viscosity of gas at temperature and pressure [Pa*s], by default None.
    heat_capacity_solid_dt_integral: Callable, optional
        A function that receives temperature1, temperature2 and pressure, and
        returns the integral of the solid heat capacity over temperature1 and
        temperature2 at pressure, by default None.
    heat_capacity_liquid_dt_integral: Callable, optional
        A function that receives temperature1, temperature2 and pressure, and
        returns the integral of the liquid heat capacity over temperature1 and
        temperature2 at pressure, by default None.
    heat_capacity_gas_dt_integral: Callable, optional
        A function that receives temperature1, temperature2 and pressure, and
        returns the integral of the gas heat capacity over temperature1 and
        temperature2 at pressure, by default None.
    vectorize_functions: bool, optional
        When True, numpy.vectorize() is applied to the temperature and pressure
        Substance Callable paramaters on Substance object init, by default
        False.

    Attributes
    ----------
    name : str, optional
       Name of the substance, by default None.
    molecular_weight : float, optional
        The molecular weight of the substance [g/mol], by default None.
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None.
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None.
    critical_temperature : float, optional
        The critical temperature of the substance [K], by default None.
    critical_pressure : float, optional
        The critical pressure of the substance [Pa], by default None.
    acentric_factor : float, optional
        The acentric factor of the substance, by default None.
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None.
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None.
    vectorize_functions: bool, optional
        When True, numpy.vectorize() is applied to the temperature and pressure
        Substance functions on Substance object init, by default False.
    """

    def __init__(
        self,
        name: str,
        molecular_weight: float = None,
        normal_boiling_point: float = None,
        normal_melting_point: float = None,
        critical_temperature: float = None,
        critical_pressure: float = None,
        acentric_factor: float = None,
        formation_enthalpy: float = None,
        formation_enthalpy_ig: float = None,
        vaporization_enthalpy: Callable = None,
        sublimation_enthalpy: Callable = None,
        volume_solid: Callable = None,
        volume_liquid: Callable = None,
        heat_capacity_solid: Callable = None,
        heat_capacity_liquid: Callable = None,
        heat_capacity_gas: Callable = None,
        thermal_conductivity_liquid: Callable = None,
        thermal_conductivity_gas: Callable = None,
        viscosity_liquid: Callable = None,
        viscosity_gas: Callable = None,
        heat_capacity_solid_dt_integral: Callable = None,
        heat_capacity_liquid_dt_integral: Callable = None,
        heat_capacity_gas_dt_integral: Callable = None,
        vectorize_functions: bool = False,
    ) -> None:
        # Symbolic init
        super().__init__(names=name)

        # Pure compound properties:
        self.name = name
        self.molecular_weight = molecular_weight
        self.normal_boiling_point = normal_boiling_point
        self.normal_melting_point = normal_melting_point
        self.critical_temperature = critical_temperature
        self.critical_pressure = critical_pressure
        self.acentric_factor = acentric_factor
        self.formation_enthalpy = formation_enthalpy
        self.formation_enthalpy_ig = formation_enthalpy_ig

        # Vectorize option
        self.vectorize_functions = vectorize_functions

        # Temperature and pressure-dependent properties calculation functions:
        # numpy vectorization is not made by default
        if self.vectorize_functions:
            self._vaporization_enthalpy = np.vectorize(
                vaporization_enthalpy, signature="()->()"
            )

            self._sublimation_enthalpy = np.vectorize(
                sublimation_enthalpy, signature="()->()"
            )

            self._volume_solid = np.vectorize(
                volume_solid, signature="(),()->()"
            )

            self._volume_liquid = np.vectorize(
                volume_liquid, signature="(),()->()"
            )

            self._heat_capacity_solid = np.vectorize(
                heat_capacity_solid, signature="(),()->()"
            )

            self._heat_capacity_liquid = np.vectorize(
                heat_capacity_liquid, signature="(),()->()"
            )

            self._heat_capacity_gas = np.vectorize(
                heat_capacity_gas, signature="(),()->()"
            )

            self._thermal_conductivity_liquid = np.vectorize(
                thermal_conductivity_liquid, signature="(),()->()"
            )

            self._thermal_conductivity_gas = np.vectorize(
                thermal_conductivity_gas, signature="(),()->()"
            )

            self._viscosity_liquid = np.vectorize(
                viscosity_liquid, signature="(),()->()"
            )

            self._viscosity_gas = np.vectorize(
                viscosity_gas, signature="(),()->()"
            )

            self._heat_capacity_solid_dt_integral = np.vectorize(
                heat_capacity_solid_dt_integral, signature="(),(),()->()"
            )

            self._heat_capacity_liquid_dt_integral = np.vectorize(
                heat_capacity_liquid_dt_integral, signature="(),(),()->()"
            )

            self._heat_capacity_gas_dt_integral = np.vectorize(
                heat_capacity_gas_dt_integral, signature="(),(),()->()"
            )

        else:
            # No numpy vectorization if the user guaranties compatibility with
            # vectors as arguments.
            self._vaporization_enthalpy = vaporization_enthalpy
            self._sublimation_enthalpy = sublimation_enthalpy
            self._volume_solid = volume_solid
            self._volume_liquid = volume_liquid
            self._heat_capacity_solid = heat_capacity_solid
            self._heat_capacity_liquid = heat_capacity_liquid
            self._heat_capacity_gas = heat_capacity_gas
            self._thermal_conductivity_liquid = thermal_conductivity_liquid
            self._thermal_conductivity_gas = thermal_conductivity_gas
            self._viscosity_liquid = viscosity_liquid
            self._viscosity_gas = viscosity_gas
            self._heat_capacity_solid_dt_integral = (
                heat_capacity_solid_dt_integral
            )
            self._heat_capacity_liquid_dt_integral = (
                heat_capacity_liquid_dt_integral
            )
            self._heat_capacity_gas_dt_integral = heat_capacity_gas_dt_integral

    @classmethod
    def from_thermo_database(
        cls, name: str, thermo_identification: str
    ) -> "Substance":
        """Substance instance from Bell Caleb's thermo library.

        Method that uses Bell Caleb's Thermo library to construct the Substance
        object.

        Cite:
        Caleb Bell and Contributors (2016-2021). Thermo: Chemical properties
        component of Chemical Engineering Design Library (ChEDL)
        https://github.com/CalebBell/thermo.

        Parameters
        ----------
        name : str
            Name that will be assigned to the Substance object.
        thermo_identification : str
            Name or CAS number of the substance that would be used to search in
            the Thermo library.

        Returns
        -------
        Substance
            Instantiated Substance object from thermo database.
        """
        return thermo_substance_constructor(cls, name, thermo_identification)

    @classmethod
    def from_pickle(cls, name_file: str) -> "Substance":
        """Read a dill Substance file and return the Substance object.

        Parameters
        ----------
        name_file : str
            Name of the file.

        Returns
        -------
        Substance : Substance
            Substance object.
        """
        with open(name_file, "rb") as f:
            return loads(f.read())

    def to_pickle(self, name_file: str) -> __file__:
        """Serialize an substance object with dill library.

        This method save an object substance as a file.

        Parameters
        ----------
        name_file : str
            Name of file to save the substance object.

        Returns
        -------
        _file_
            A binary file with substance predefine object.
        """
        with open(name_file, "wb") as f:
            f.write(dumps(self))

    def vaporization_enthalpy(self, temperature: float) -> float:
        """Return the vaporization enthalpy at a given temperature.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]

        Returns
        -------
        float or ndarray[float]
            Vaporization enthalpy. [J/mol]
        """
        return self._vaporization_enthalpy(temperature)

    def sublimation_enthalpy(self, temperature: float) -> float:
        """Return the sublimation enthalpy at a given temperature.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]

        Returns
        -------
        float or ndarray[float]
            Sublimation enthalpy. [J/mol]
        """
        return self._sublimation_enthalpy(temperature)

    def fusion_enthalpy(self, temperature: float) -> float:
        """Return the fusion enthalpy at a given temperature.

        Uses the sublimation and vaporization enthalpy functions for the fusion
        enthalpy calculations at a given temperature, by calculating the
        sublimation and vaporization enthalpy difference.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]

        Returns
        -------
        float or ndarray[float]
            Fusion enthalpy. [J/mol]
        """
        fusion_h = self._sublimation_enthalpy(
            temperature
        ) - self._vaporization_enthalpy(temperature)
        return fusion_h

    def volume_solid(self, temperature: float, pressure: float) -> float:
        """Return the solid molar volume at a given temperature and pressure.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Solid molar volume. [m³/mol]
        """
        return self._volume_solid(temperature, pressure)

    def volume_liquid(self, temperature: float, pressure: float) -> float:
        """Return the liquid molar volume at a given temperature and pressure.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Liquid molar volume. [m³/mol]
        """
        return self._volume_liquid(temperature, pressure)

    def heat_capacity_solid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the pure solid heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Solid heat capacity. [J/mol/K]
        """
        return self._heat_capacity_solid(temperature, pressure)

    def heat_capacity_liquid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the pure liquid heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Liquid heat capacity. [J/mol/K]
        """
        return self._heat_capacity_liquid(temperature, pressure)

    def heat_capacity_gas(self, temperature: float, pressure: float) -> float:
        """Return the ideal gas heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Ideal gas heat capacity. [J/mol/K]
        """
        return self._heat_capacity_gas(temperature, pressure)

    def thermal_conductivity_liquid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the liquid thermal conductivity.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Liquid thermal conductivity. [W/m/K]
        """
        return self._thermal_conductivity_liquid(temperature, pressure)

    def thermal_conductivity_gas(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the gas thermal conductivity.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Gas thermal conductivity. [W/m/K]
        """
        return self._thermal_conductivity_gas(temperature, pressure)

    def viscosity_liquid(self, temperature: float, pressure: float) -> float:
        """Return the pure liquid viscosity.

        Pure liquid viscosity at a given temperature and pressure.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Liquid viscosity. [Pa s]
        """
        return self._viscosity_liquid(temperature, pressure)

    def viscosity_gas(self, temperature: float, pressure: float) -> float:
        """Return the pure gas viscosity.

         Pure gas viscosity at a given temperature and pressure.

        Parameters
        ----------
        temperature : float or ndarray[float]
            Temperature in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Gas viscosity. [Pa s]
        """
        return self._viscosity_gas(temperature, pressure)

    def heat_capacity_solid_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of solid heat capacity.

        Calculate the integral of solid heat capacity between temperature1
        and temperature2

        .. math::
            \int_{T_1}^{T_2} {C_{p_{solid}} (T, P)} \mathrm{d}T

        | :math:`T_1`: temperature1.
        | :math:`T_2`: temperature2.
        | :math:`C_{p_{solid}} (T, P)`: solid heat capacity.

        Parameters
        ----------
        temperature1 : float or ndarray[float]
            Lower temperature integral bound in Kelvin degrees. [K]
        temperature2 : float or ndarray[float]
            Upper temperature integral bound in Kelvin degrees. [K]
        pressure: float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Integral of solid heat capacity between temperature1 and
            temperature2. [J/mol]
        """
        integral = self._heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
        return integral

    def heat_capacity_liquid_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of liquid heat capacity.

        Calculate the definite integral of liquid heat capacity between
        temperature1 and temperature2

        .. math::
            \int_{T_1}^{T_2} {C_{p_{liquid}} (T, P)} \mathrm{d}T

        | :math:`T_1`: temperature1.
        | :math:`T_2`: temperature2.
        | :math:`C_{p_{liquid}} (T, P)`: liquid heat capacity.

        Parameters
        ----------
        temperature1 : float or ndarray[float]
            Lower temperature integral bound in Kelvin degrees [K]
        temperature2 : float or ndarray[float]
            Upper temperature integral bound in Kelvin degrees [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]heat_capacity_gas_dt_integral
            Definite integral of liquid heat capacity between temperature1
            and temperature2. [J/mol]
        """
        integral = self._heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
        return integral

    def heat_capacity_gas_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of gas heat capacity.

        Calculate the definite integral of gas heat capacity between
        temperature1 and temperature2.

        .. math::
            \int_{T_1}^{T_2} {C_{p_{gas}} (T, P)} \mathrm{d}T

        | :math:`T_1`: temperature1.
        | :math:`T_2`: temperature2.
        | :math:`C_{p_{gas}} (T, P)`: ideal gas heat capacity.

        Parameters
        ----------
        temperature1 : float or ndarray[float]
            Lower temperature integral bound in Kelvin degrees. [K]
        temperature2 : float or ndarray[float]
            Upper temperature integral bound in Kelvin degrees. [K]
        pressure : float or ndarray[float]
            Pressure in Pascal. [Pa]

        Returns
        -------
        float or ndarray[float]
            Definite integral of gas heat capacity between temperature1 and
            temperature2. [J/mol]
        """
        integral = self._heat_capacity_gas_dt_integral(
            temperature1, temperature2, pressure
        )
        return integral
