"""Substance module.

Class to define a substance for ReactorD library.
"""
from typing import Callable

from scipy.integrate import quad

from thermo.chemical import Chemical


class Substance:
    """Substance object class.

    Class to define a substance object. Specific attributes definition will be
    required for the reactors, described in each reactor documentation. For
    example, an adiabatic reactor will require that substances define a heat
    capacity function, on the other hand, when using isothermic reactors will
    not be necessary. Substance has the .from_thermo_data_base alternative
    construction method.
    E.g:

    water = Substance.from_thermo_data_base('water')

    Parameters
    ----------
    name : str, optional
       Name of the substance, by default None
    mw : float, optional
        The molecular weight of the substance [g/mol], by default None
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None
    tc : float, optional
        The critical temperature of the substance [K], by default None
    pc : float, optional
        The critical pressure of the substance [Pa], by default None
    omega : float, optional
        The acentric factor of the substance, by default None
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    formation_gibbs : float, optional
       Standard state molar change of Gibbs energy of formation [J/mol], by
       default None
    formation_gibbs_ig : float, optional
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by default
        None
    vaporization_enthalpy_t : Callable, optional
        A function that receives temperature and returns the vaporization
        enthalpy at temperature [J/mol], by default None
    sublimation_enthalpy_t : Callable, optional
        A function that receives temperature and returns the sublimation
        enthalpy at temperature [J/mol], by default None
    volume_s_t : Callable, optional
        A function that receives temperature and returns the molar volume of
        the solid at temperature [m^3/mol], by default None
    volume_l_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        molar volume of liquid at temperature and pressure [m^3/mol], by
        default None
    volume_g_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        molar volume of gas at temperature and pressure [m^3/mol], by default
        None
    heat_capacity_s_t : Callable, optional
        A function that receives temperature and pressure, and returns the heat
        capacity of solid at temperature and pressure [J/mol/K], by default
        None
    heat_capacity_l_t : Callable, optional
        A function that receives temperature and returns the heat capacity of
        liquid at temperature [J/mol/K], by default None
    heat_capacity_g_t : Callable, optional
        A function that receives temperature and returns the heat capacity of
        gas at temperature [J/mol/K], by default None
    thermal_conductivity_l_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        thermal conductivity of liquid at temperature and pressure [W/m/K], by
        default None
    thermal_conductivity_g_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        thermal conductivity of gas at temperature and pressure [W/m/K], by
        default None
    viscosity_l_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        viscosity of liquid at temperature and pressure [Pa*s], by default None
    viscosity_g_tp : Callable, optional
        A function that receives temperature and pressure, and returns the
        viscosity of gas at temperature and pressure [Pa*s], by default None

    Attributes
    ----------
    name : str, optional
       Name of the substance, by default None
    mw : float, optional
        The molecular weight of the substance [g/mol], by default None
    normal_boiling_point : float, optional
        The normal boiling point of the substance [K], by default None
    normal_melting_point : float, optional
        The normal melting point of the substance [K], by default None
    tc : float, optional
        The critical temperature of the substance [K], by default None
    pc : float, optional
        The critical pressure of the substance [Pa], by default None
    omega : float, optional
        The acentric factor of the substance, by default None
    formation_enthalpy : float, optional
        Standard state molar enthalpy of formation [J/mol], by default None
    formation_enthalpy_ig : float, optional
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    formation_gibbs : float, optional
       Standard state molar change of Gibbs energy of formation [J/mol], by
       default None
    formation_gibbs_ig : float, optional
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by default
        None
    """

    def __init__(
        self,
        name: str = None,
        mw: float = None,
        normal_boiling_point: float = None,
        normal_melting_point: float = None,
        tc: float = None,
        pc: float = None,
        omega: float = None,
        formation_enthalpy: float = None,
        formation_enthalpy_ig: float = None,
        formation_gibbs: float = None,
        formation_gibbs_ig: float = None,
        vaporization_enthalpy_t: Callable = None,
        sublimation_enthalpy_t: Callable = None,
        volume_s_t: Callable = None,
        volume_l_tp: Callable = None,
        volume_g_tp: Callable = None,
        heat_capacity_s_t: Callable = None,
        heat_capacity_l_t: Callable = None,
        heat_capacity_g_t: Callable = None,
        thermal_conductivity_l_tp: Callable = None,
        thermal_conductivity_g_tp: Callable = None,
        viscosity_l_tp: Callable = None,
        viscosity_g_tp: Callable = None,
    ) -> None:

        # Pure compound properties:
        self.name = name
        self.mw = mw
        self.normal_boiling_point = normal_boiling_point
        self.normal_melting_point = normal_melting_point
        self.tc = tc
        self.pc = pc
        self.omega = omega
        self.formation_enthalpy = formation_enthalpy
        self.formation_enthalpy_ig = formation_enthalpy_ig
        self.formation_gibbs = formation_gibbs
        self.formation_gibbs_ig = formation_gibbs_ig
        # Temperature-dependent properties calculation functions:
        self._vaporization_enthalpy_t = vaporization_enthalpy_t
        self._sublimation_enthalpy_t = sublimation_enthalpy_t
        self._volume_s_t = volume_s_t
        self._volume_l_tp = volume_l_tp
        self._volume_g_tp = volume_g_tp
        self._heat_capacity_s_t = heat_capacity_s_t
        self._heat_capacity_l_t = heat_capacity_l_t
        self._heat_capacity_g_t = heat_capacity_g_t
        # self._thermal_conductivity_s_tp = thermal_conductivity_s_tp
        self._thermal_conductivity_l_tp = thermal_conductivity_l_tp
        self._thermal_conductivity_g_tp = thermal_conductivity_g_tp
        self._viscosity_l_tp = viscosity_l_tp
        self._viscosity_g_tp = viscosity_g_tp

    @classmethod
    def from_thermo_database(cls, identification: str):
        """Substance instance from Bell Caleb's thermo library.

        Method that uses Bell Caleb's Thermo library to construct the Substance
        object.

        Cite:

        Caleb Bell and Contributors (2016-2021). Thermo: Chemical properties
        component of Chemical Engineering Design Library (ChEDL)
        https://github.com/CalebBell/thermo.

        Parameters
        ----------
        identification : string
            Name or CAS number of the substance

        Returns
        -------
        Substance
            Instantiated Substance object from thermo database.
        """
        chemobj = Chemical(identification)

        substance_object = cls(
            name=chemobj.name,
            mw=chemobj.MW,
            normal_boiling_point=chemobj.Tb,
            normal_melting_point=chemobj.Tm,
            tc=chemobj.Tc,
            pc=chemobj.Pc,
            omega=chemobj.omega,
            formation_enthalpy=chemobj.Hfm,
            formation_enthalpy_ig=chemobj.Hfgm,
            formation_gibbs=chemobj.Gfm,
            formation_gibbs_ig=chemobj.Gfgm,
            vaporization_enthalpy_t=chemobj.EnthalpyVaporization,
            sublimation_enthalpy_t=chemobj.EnthalpySublimation,
            volume_s_t=chemobj.VolumeSolid,
            volume_l_tp=chemobj.VolumeLiquid,
            volume_g_tp=chemobj.VolumeGas,
            heat_capacity_s_t=chemobj.HeatCapacitySolid,
            heat_capacity_l_t=chemobj.HeatCapacityLiquid,
            heat_capacity_g_t=chemobj.HeatCapacityGas,
            thermal_conductivity_l_tp=chemobj.ThermalConductivityLiquid,
            thermal_conductivity_g_tp=chemobj.ThermalConductivityGas,
            viscosity_l_tp=chemobj.ViscosityLiquid,
            viscosity_g_tp=chemobj.ViscosityGas,
        )
        return substance_object

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
        return self._vaporization_enthalpy_t(temperature)

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
        return self._sublimation_enthalpy_t(temperature)

    def fusion_enthalpy(self, temperature: float) -> float:
        """Return the fusion enthalpy at a given temperature.

        Uses the sublimation and vaporization enthalpy functions for the fusion
        enthalpy calculations at a given temperature, by calculating the
        sublimation and vaporization enthalpy difference.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Fusion enthalpy in Joule per mol [J/mol]
        """
        fusion_h = self._sublimation_enthalpy_t(
            temperature
        ) - self._vaporization_enthalpy_t(temperature)
        return fusion_h

    def volume_solid(self, temperature: float) -> float:
        """Return the solid molar volume at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Solid molar volume in cubic meters per mol [m^3/mol]
        """
        return self._volume_s_t(temperature)

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
            Liquid molar volume in cubic meters per mol [m^3/mol]
        """
        return self._volume_l_tp(temperature, pressure)

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
            Liquid molar volume in cubic meters per mol [m^3/mol]
        """
        return self._volume_g_tp(temperature, pressure)

    def heat_capacity_solid(self, temperature: float) -> float:
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
        return self._heat_capacity_s_t(temperature)

    def heat_capacity_liquid(self, temperature: float) -> float:
        """Return the pure liquid heat capacity at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin degrees [K]

        Returns
        -------
        float
            Pure liquid heat capacity in Joule per mol per Kelvin [J/mol/K]
        """
        return self._heat_capacity_l_t(temperature)

    def heat_capacity_gas(self, temperature: float) -> float:
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
        return self._heat_capacity_g_t(temperature)

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
            Liquid thermal conductivity in Watts per meter per Kelvin [W/m/K]
        """
        return self._thermal_conductivity_l_tp(temperature, pressure)

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
            Gas thermal conductivity in Watts per meter per Kelvin [W/m/K]
        """
        return self._thermal_conductivity_g_tp(temperature, pressure)

    def viscosity_liquid(self, temperature: float, pressure: float) -> float:
        """Return the pure liquid viscosity.

         Pure liquid viscosity as a function of temperature and pressure

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
        return self._viscosity_l_tp(temperature, pressure)

    def viscosity_gas(self, temperature: float, pressure: float) -> float:
        """Return the pure gas viscosity.

         Pure gas viscosity as a function of temperature and pressure

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
        return self._viscosity_g_tp(temperature, pressure)

    def heat_capacity_solid_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:
        """Return the integral of solid heat capacity.

        Calculate the integral of solid heat capacity between temperature1 and
        temperature2

        Parameters
        ----------
        temperature1 : float
            Temperature 1 in Kelvin degrees [K]
        temperature2 : float
            Temperature 2 in Kelvin degrees [K]

        Returns
        -------
        float
            Integral of solid heat capacity between temperature1 and
            temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_solid, a=temperature1, b=temperature2
        )

        return integral

    def heat_capacity_liquid_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:
        """Return the integral of liquid heat capacity.

        Calculate the definite integral of liquid heat capacity between
        temperature1 and temperature2

        Parameters
        ----------
        temperature1 : float
            Temperature 1 in Kelvin degrees [K]
        temperature2 : float
            Temperature 2 in Kelvin degrees [K]

        Returns
        -------
        float
            Definite integral of liquid heat capacity between temperature1
            and temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_liquid, a=temperature1, b=temperature2
        )

        return integral

    def heat_capacity_gas_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:
        """Return the integral of gas heat capacity.

        Calculate the definite integral of gas heat capacity between
        temperature1 and temperature2

        Parameters
        ----------
        temperature1 : float
            Temperature 1 in Kelvin degrees [K]
        temperature2 : float
            Temperature 2 in Kelvin degrees [K]

        Returns
        -------
        float
            Definite integral of gas heat capacity between temperature1
            and temperature2 in Joule per mol per Kelvin [J/mol/K]
        """
        integral, err = quad(
            self.heat_capacity_gas, a=temperature1, b=temperature2
        )

        return integral
