"""Substance module.

Class to define a substance for ReactorD library.
"""
from typing import Callable

from dill import dumps, loads

import numpy as np

from scipy.integrate import quad

from thermo import ChemicalConstantsPackage


class Substance:
    """Substance object class.

    Class to define a substance object. Specific attributes definition will be
    required for the reactors, described in each reactor documentation. For
    example, an adiabatic reactor will require that substances define a heat
    capacity function, on the other hand, when using isothermic reactors this
    won't be necessary. Substance has the .from_thermo_data_base alternative
    construction method.
    E.g:

    water = Substance.from_thermo_data_base('water')

    Parameters
    ----------
    name : str, optional
       Name of the substance, by default None
    molecular_weight : float, optional
        The molecular weight of the substance [g/mol], by default None
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None
    critical_temperature : float, optional
        The critical temperature of the substance [K], by default None
    critical_pressure : float, optional
        The critical pressure of the substance [Pa], by default None
    acentric_factor : float, optional
        The acentric factor of the substance, by default None
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    formation_gibbs_ig : float, optional
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by default
        None
    vaporization_enthalpy : Callable, optional
        A function that receives a temperature and returns the vaporization
        enthalpy at that temperature [J/mol], by default None
    sublimation_enthalpy : Callable, optional
        A function that receives a temperature and returns the sublimation
        enthalpy at that temperature [J/mol], by default None
    volume_solid : Callable, optional
        A function that receives a temperature and returns the molar volume of
        the solid at that temperature [m³/mol], by default None
    volume_liquid : Callable, optional
        A function that receives a temperature and pressure, and returns the
        molar volume of liquid at that temperature and pressure [m³/mol], by
        default None
    volume_gas : Callable, optional
        A function that receives a temperature and pressure, and returns the
        molar volume of the gas at that temperature and pressure [m³/mol],
        by default None
    heat_capacity_solid : Callable, optional
        A function that receives a temperature and pressure, and returns the
        heat capacity of the solid at that temperature and pressure [J/mol/K],
        by default None
    heat_capacity_liquid : Callable, optional
        A function that receives a temperature and returns the heat capacity of
        the liquid at that temperature [J/mol/K], by default None
    heat_capacity_gas : Callable, optional
        A function that receives a temperature and returns the heat capacity of
        the gas at that temperature [J/mol/K], by default None
    thermal_conductivity_liquid : Callable, optional
        A function that receives a temperature and pressure, and returns the
        thermal conductivity of the liquid at that temperature and pressure
        [W/m/K], by default None
    thermal_conductivity_gas : Callable, optional
        A function that receives a temperature and pressure, and returns the
        thermal conductivity of the gas at that temperature and pressure
        [W/m/K], by default None
    viscosity_liquid : Callable, optional
        A function that receives a temperature and pressure, and returns the
        viscosity of thr liquid at that temperature and pressure [Pa*s], by
        default None
    viscosity_gas : Callable, optional
        A function that receives temperature and pressure, and returns the
        viscosity of gas at temperature and pressure [Pa*s], by default None

    Attributes
    ----------
    name : str, optional
       Name of the substance, by default None
    molecular_weight : float, optional
        The molecular weight of the substance [g/mol], by default None
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None
    critical_temperature : float, optional
        The critical temperature of the substance [K], by default None
    critical_pressure : float, optional
        The critical pressure of the substance [Pa], by default None
    acentric_factor : float, optional
        The acentric factor of the substance, by default None
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    formation_gibbs_ig : float, optional
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by default
        None
    """

    def __init__(
        self,
        name: str = None,
        molecular_weight: float = None,
        normal_boiling_point: float = None,
        normal_melting_point: float = None,
        critical_temperature: float = None,
        critical_pressure: float = None,
        acentric_factor: float = None,
        formation_enthalpy: float = None,
        formation_enthalpy_ig: float = None,
        formation_gibbs_ig: float = None,
        vaporization_enthalpy: Callable = None,
        sublimation_enthalpy: Callable = None,
        volume_solid: Callable = None,
        volume_liquid: Callable = None,
        volume_gas: Callable = None,
        heat_capacity_solid: Callable = None,
        heat_capacity_liquid: Callable = None,
        heat_capacity_gas: Callable = None,
        thermal_conductivity_liquid: Callable = None,
        thermal_conductivity_gas: Callable = None,
        viscosity_liquid: Callable = None,
        viscosity_gas: Callable = None,
    ) -> None:

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
        self.formation_gibbs_ig = formation_gibbs_ig

        # Temperature-dependent properties calculation functions:
        self._vaporization_enthalpy = np.vectorize(
            vaporization_enthalpy, signature="()->()"
        )
        self._sublimation_enthalpy = np.vectorize(
            sublimation_enthalpy, signature="()->()"
        )
        self._volume_solid = np.vectorize(volume_solid, signature="(),()->()")
        self._volume_liquid = np.vectorize(
            volume_liquid, signature="(),()->()"
        )
        self._volume_gas = np.vectorize(volume_gas, signature="(),()->()")
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

    @classmethod
    def from_thermo_database(cls, identification: str):
        """Substance instance from Bell Caleb's thermo library.

        Method that uses Bell Caleb's Thermo library to construct the
        Substance object.

        Cite:

        Caleb Bell and Contributors (2016-2021). Thermo: Chemical
        properties component of Chemical Engineering Design
        Library (ChEDL) https://github.com/CalebBell/thermo.

        Parameters
        ----------
        identification : string
            Name or CAS number of the substance

        Returns
        -------
        Substance
            Instantiated Substance object from thermo database.
        """
        corr = ChemicalConstantsPackage.correlations_from_IDs([identification])

        # Temperature functions
        def vaporization_enthalpy(temperature: float) -> float:
            enthalpy = corr.EnthalpyVaporizations[0].T_dependent_property(
                temperature
            )
            return enthalpy

        def sublimation_enthalpy(temperature: float) -> float:
            enthalpy = corr.EnthalpySublimations[0].T_dependent_property(
                temperature
            )
            return enthalpy

        # Temperature and pressure dependent functions
        def volume_solid(temperature: float, pressure: float) -> float:
            volume = corr.VolumeSolids[0].T_dependent_property(temperature)
            return volume

        def volume_liquid(temperature: float, pressure: float) -> float:
            volume = corr.VolumeLiquids[0].T_dependent_property(temperature)
            return volume

        def volume_gas(temperature: float, pressure: float) -> float:
            method = corr.VolumeGases[0].method_P
            volume = corr.VolumeGases[0].calculate_P(
                temperature, pressure, method
            )
            return volume

        def heat_capacity_solid(temperature: float, pressure: float) -> float:
            heat_cap = corr.HeatCapacitySolids[0].T_dependent_property(
                temperature
            )
            return heat_cap

        def heat_capacity_liquid(temperature: float, pressure: float) -> float:
            heat_cap = corr.HeatCapacityLiquids[0].T_dependent_property(
                temperature
            )
            return heat_cap

        def heat_capacity_gas(temperature: float, pressure: float) -> float:
            heat_cap = corr.HeatCapacityGases[0].T_dependent_property(
                temperature
            )
            return heat_cap

        def thermal_conductivity_liquid(
            temperature: float, pressure: float
        ) -> float:
            thermal_cond = corr.ThermalConductivityLiquids[
                0
            ].T_dependent_property(temperature)
            return thermal_cond

        def thermal_conductivity_gas(
            temperature: float, pressure: float
        ) -> float:
            thermal_cond = corr.ThermalConductivityGases[
                0
            ].T_dependent_property(temperature)
            return thermal_cond

        def viscosity_liquid(temperature: float, pressure: float) -> float:
            thermal_cond = corr.ViscosityLiquids[0].T_dependent_property(
                temperature
            )
            return thermal_cond

        def viscosity_gas(temperature: float, pressure: float) -> float:
            thermal_cond = corr.ViscosityGases[0].T_dependent_property(
                temperature
            )
            return thermal_cond

        substance_object = cls(
            name=corr.constants.names[0],
            molecular_weight=corr.constants.MWs[0],
            normal_boiling_point=corr.constants.Tbs[0],
            normal_melting_point=corr.constants.Tms[0],
            critical_temperature=corr.constants.Tcs[0],
            critical_pressure=corr.constants.Pcs[0],
            acentric_factor=corr.constants.omegas[0],
            formation_enthalpy=corr.constants.Hf_STPs[0],
            formation_enthalpy_ig=corr.constants.Hfgs[0],
            formation_gibbs_ig=corr.constants.Gfgs[0],
            vaporization_enthalpy=vaporization_enthalpy,
            sublimation_enthalpy=sublimation_enthalpy,
            volume_solid=volume_solid,
            volume_liquid=volume_liquid,
            volume_gas=volume_gas,
            heat_capacity_solid=heat_capacity_solid,
            heat_capacity_liquid=heat_capacity_liquid,
            heat_capacity_gas=heat_capacity_gas,
            thermal_conductivity_liquid=thermal_conductivity_liquid,
            thermal_conductivity_gas=thermal_conductivity_gas,
            viscosity_liquid=viscosity_liquid,
            viscosity_gas=viscosity_gas,
        )
        return substance_object

    @classmethod
    def from_pickle(cls, name_file: str) -> "Substance":
        """Read a dill Substance file and return the Substance object.

        Cite:
        https://github.com/uqfoundation/dill

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

        Cite:
        https://github.com/uqfoundation/dill

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
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Vaporization enthalpy in Joule per mol [J/mol]
        """
        return self._vaporization_enthalpy(temperature)

    def sublimation_enthalpy(self, temperature: float) -> float:
        """Return the sublimation enthalpy at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Sublimation enthalpy in Joule per mol [J/mol]
        """
        return self._sublimation_enthalpy(temperature)

    def fusion_enthalpy(self, temperature: float) -> float:
        """Return the fusion enthalpy at a given temperature.

        Uses the sublimation and vaporization enthalpy functions for the
        fusion enthalpy calculations at a given temperature, by
        calculating the sublimation and vaporization enthalpy difference.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Fusion enthalpy in Joule per mol [J/mol]
        """
        fusion_h = self._sublimation_enthalpy(
            temperature
        ) - self._vaporization_enthalpy(temperature)
        return fusion_h

    def volume_solid(self, temperature: float, pressure: float) -> float:
        """Return the solid molar volume at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Solid molar volume in cubic meters per mol [m³/mol]
        """
        return self._volume_solid(temperature, pressure)

    def volume_liquid(self, temperature: float, pressure: float) -> float:
        """Return the liquid molar volume at a given temperature and pressure.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in pascal [Pa]

        Returns
        -------
        float
            Liquid molar volume in cubic meters per mol [m³/mol]
        """
        return self._volume_liquid(temperature, pressure)

    def volume_gas(self, temperature: float, pressure: float) -> float:
        """Return the gas molar volume at a given temperature and pressure.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in Pascal [Pa]

        Returns
        -------
        float
            Liquid molar volume in cubic meters per mol [m³/mol]
        """
        return self._volume_gas(temperature, pressure)

    def heat_capacity_solid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the pure solid heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Pure solid heat capacity in Joule per mol per Kelvin [J/mol/K]
        """
        return self._heat_capacity_solid(temperature, pressure)

    def heat_capacity_liquid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the pure liquid heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Pure liquid heat capacity in Joule per mol per Kelvin
            [J/mol/K]
        """
        return self._heat_capacity_liquid(temperature, pressure)

    def heat_capacity_gas(self, temperature: float, pressure: float) -> float:
        """Return the pure gas heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Pure gas heat capacity in Joule per mol per Kelvin [J/mol/K]
        """
        return self._heat_capacity_gas(temperature, pressure)

    def thermal_conductivity_liquid(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the liquid thermal conductivity.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in Pascal [Pa]

        Returns
        -------
        float
            Liquid thermal conductivity in Watts per meter per Kelvin
            [W/m/K]
        """
        return self._thermal_conductivity_liquid(temperature, pressure)

    def thermal_conductivity_gas(
        self, temperature: float, pressure: float
    ) -> float:
        """Return the gas thermal conductivity.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in Pascal [Pa]

        Returns
        -------
        float
            Gas thermal conductivity in Watts per meter per Kelvin
            [W/m/K]
        """
        return self._thermal_conductivity_gas(temperature, pressure)

    def viscosity_liquid(self, temperature: float, pressure: float) -> float:
        """Return the pure liquid viscosity.

        Pure liquid viscosity at a given temperature and pressure.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in Pascal [Pa]

        Returns
        -------
        float
            Pure liquid viscosity in [Pa*s]
        """
        return self._viscosity_liquid(temperature, pressure)

    def viscosity_gas(self, temperature: float, pressure: float) -> float:
        """Return the pure gas viscosity.

         Pure gas viscosity at a given temperature and pressure.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]
        pressure : float
            Pressure in Pascal [Pa]

        Returns
        -------
        float
            Pure gas viscosity in [Pa*s]
        """
        return self._viscosity_gas(temperature, pressure)

    def heat_capacity_solid_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of solid heat capacity.

        Calculate the integral of solid heat capacity between temperature1
        and temperature2

        .. math::
            \int_{T_1}^{T_2} {C_{p_{solid}} (T)} \mathrm{d}T

        Parameters
        ----------
        temperature1 : float
            Lower temperature integral bound in Kelvin degrees [K]
        temperature2 : float
            Upper temperature integral bound in Kelvin degrees [K]

        Returns
        -------
        float
            Integral of solid heat capacity between temperature1 and
            temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_solid,
            a=temperature1,
            b=temperature2,
            args=(pressure,),
        )

        return integral

    def heat_capacity_liquid_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of liquid heat capacity.

        Calculate the definite integral of liquid heat capacity between
        temperature1 and temperature2

        .. math::
            \int_{T_1}^{T_2} {C_{p_{liquid}} (T)} \mathrm{d}T

        Parameters
        ----------
        temperature1 : float
            Lower temperature integral bound in Kelvin degrees [K]
        temperature2 : float
            Upper temperature integral bound in Kelvin degrees [K]

        Returns
        -------
        float
            Definite integral of liquid heat capacity between temperature1
            and temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_liquid,
            a=temperature1,
            b=temperature2,
            args=(pressure,),
        )

        return integral

    def heat_capacity_gas_dt_integral(
        self, temperature1: float, temperature2: float, pressure: float
    ) -> float:
        r"""Return the integral of gas heat capacity.

        Calculate the definite integral of gas heat capacity between
        temperature1 and temperature2

        .. math::
            \int_{T_1}^{T_2} {C_{p_{gas}} (T)} \mathrm{d}T

        Parameters
        ----------
        temperature1 : float
            Lower temperature integral bound in Kelvin degrees [K]
        temperature2 : float
            Upper temperature integral bound in Kelvin degrees [K]

        Returns
        -------
        float
            Definite integral of gas heat capacity between temperature1
            and temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_gas,
            a=temperature1,
            b=temperature2,
            args=(pressure,),
        )

        return integral
