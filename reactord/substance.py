"""Substance module.

Class to define a substance for ReactorD library.
"""
import pickle
from typing import Callable

from reactord.utils import vectorize

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
    formation_gibbs : float, optional
       Standard state molar change of Gibbs energy of formation [J/mol], by
       default None
    formation_gibbs_ig : float, optional
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by default
        None
    vaporization_enthalpy_t : Callable, optional
        A function that receives a temperature and returns the vaporization
        enthalpy at that temperature [J/mol], by default None
    sublimation_enthalpy_t : Callable, optional
        A function that receives a temperature and returns the sublimation
        enthalpy at that temperature [J/mol], by default None
    volume_solid_t : Callable, optional
        A function that receives a temperature and returns the molar volume of
        the solid at that temperature [m³/mol], by default None
    volume_liquid_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        molar volume of liquid at that temperature and pressure [m³/mol], by
        default None
    volume_gas_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        molar volume of the gas at that temperature and pressure [m^3/mol],
        by default None
    heat_capacity_solid_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        heat capacity of the solid at that temperature and pressure [J/mol/K],
        by default None
    heat_capacity_liquid_t : Callable, optional
        A function that receives a temperature and returns the heat capacity of
        the liquid at that temperature [J/mol/K], by default None
    heat_capacity_gas_t : Callable, optional
        A function that receives a temperature and returns the heat capacity of
        the gas at that temperature [J/mol/K], by default None
    thermal_conductivity_liquid_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        thermal conductivity of the liquid at that temperature and pressure
        [W/m/K], by default None
    thermal_conductivity_gas_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        thermal conductivity of the gas at that temperature and pressure
        [W/m/K], by default None
    viscosity_liquid_tp : Callable, optional
        A function that receives a temperature and pressure, and returns the
        viscosity of thr liquid at that temperature and pressure [Pa*s], by
        default None
    viscosity_gas_tp : Callable, optional
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
        molecular_weight: float = None,
        normal_boiling_point: float = None,
        normal_melting_point: float = None,
        critical_temperature: float = None,
        critical_pressure: float = None,
        acentric_factor: float = None,
        formation_enthalpy: float = None,
        formation_enthalpy_ig: float = None,
        formation_gibbs: float = None,
        formation_gibbs_ig: float = None,
        vaporization_enthalpy_t: Callable = None,
        sublimation_enthalpy_t: Callable = None,
        volume_solid_t: Callable = None,
        volume_liquid_tp: Callable = None,
        volume_gas_tp: Callable = None,
        heat_capacity_solid_tp: Callable = None,
        heat_capacity_liquid_t: Callable = None,
        heat_capacity_gas_t: Callable = None,
        thermal_conductivity_liquid_tp: Callable = None,
        thermal_conductivity_gas_tp: Callable = None,
        viscosity_liquid_tp: Callable = None,
        viscosity_gas_tp: Callable = None,
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
        self.formation_gibbs = formation_gibbs
        self.formation_gibbs_ig = formation_gibbs_ig
        # Temperature-dependent properties calculation functions:
        self._vaporization_enthalpy_t = vaporization_enthalpy_t
        self._sublimation_enthalpy_t = sublimation_enthalpy_t
        self._volume_solid_t = volume_solid_t
        self._volume_liquid_tp = volume_liquid_tp
        self._volume_gas_tp = volume_gas_tp
        self._heat_capacity_solid_tp = heat_capacity_solid_tp
        self._heat_capacity_liquid_t = heat_capacity_liquid_t
        self._heat_capacity_gas_t = heat_capacity_gas_t
        # self._thermal_conductivity_s_tp = thermal_conductivity_s_tp
        self._thermal_conductivity_liquid_tp = thermal_conductivity_liquid_tp
        self._thermal_conductivity_gas_tp = thermal_conductivity_gas_tp
        self._viscosity_liquid_tp = viscosity_liquid_tp
        self._viscosity_gas_tp = viscosity_gas_tp

    def create_substance_file(self, name_file) -> __file__:
        """Serialize an object substance.

        This method save an object substance as a file

        Parameters
        ----------
        name_file : str
            name of file to save the substance object

        Returns
        -------
        _file_
            A binary file with substance predefine object
        """
        with open(name_file, "wb") as f:
            return pickle.dump(self, f)

    @classmethod
    def load_file(cls, name_file):
        """Read save file of substance object.

        Parameters
        ----------
        name_file : str

        Returns
        -------
        Substance : Substance
        """
        with open(name_file, "rb") as f:
            return pickle.load(f)

    @classmethod
    def from_thermo_database(cls, identification: str):
        """Substance instance from Bell Caleb's thermo library.

        Method that uses Bell Caleb's Thermo library to construct the
        Substance object.

        Cite:

        Caleb Bell and Contributors (2016-2021). Thermo: Chemical
        properties        component of Chemical Engineering Design
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

        # pressure dependent functions



        substance_object = cls(
            name=corr.constants.names[0],
            molecular_weight=corr.constants.Mws[0],
            normal_boiling_point=corr.constants.Tbs[0],
            normal_melting_point=corr.constants.Tms[0],
            critical_temperature=corr.constants.Tcs[0],
            critical_pressure=corr.constants.Pcs[0],
            acentric_factor=corr.constants.omegas[0],
            formation_enthalpy=corr.constants.Hf_STPs[0],
            formation_enthalpy_ig=corr.constants.Hfgs[0],
            formation_gibbs=0, # TODO look Hvap_298s
            formation_gibbs_ig=corr.constants.Gfgs[0],
            vaporization_enthalpy_t=corr.EnthalpyVaporization[0].T_dependent_property,
            sublimation_enthalpy_t=corr.EnthalpySublimation,
            volume_solid_tp=corr.VolumeSolid,
            volume_liquid_tp=corr.VolumeLiquid,
            volume_gas_tp=corr.VolumeGas,
            heat_capacity_solid_tp=corr.HeatCapacitySolid,
            heat_capacity_liquid_tp=corr.HeatCapacityLiquid,
            heat_capacity_gas_tp=corr.HeatCapacityGas,
            thermal_conductivity_liquid_tp=corr.ThermalConductivityLiquid,
            thermal_conductivity_gas_tp=corr.ThermalConductivityGas,
            viscosity_liquid_tp=corr.ViscosityLiquid,
            viscosity_gas_tp=corr.ViscosityGas,
        )
        return substance_object

    @vectorize(signature="()->()", excluded={0})
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

    @vectorize(signature="()->()", excluded={0})
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

    @vectorize(signature="()->()", excluded={0})
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
        fusion_h = self._sublimation_enthalpy_t(
            temperature
        ) - self._vaporization_enthalpy_t(temperature)
        return fusion_h

    @vectorize(signature="()->()", excluded={0})
    def volume_solid(self, temperature: float) -> float:
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
        return self._volume_solid_t(temperature)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._volume_liquid_tp(temperature, pressure)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._volume_gas_tp(temperature, pressure)

    @vectorize(signature="()->()", excluded={0})
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
        return self._heat_capacity_solid_tp(temperature)

    @vectorize(signature="()->()", excluded={0})
    def heat_capacity_liquid(self, temperature: float) -> float:
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
        return self._heat_capacity_liquid_t(temperature)

    @vectorize(signature="()->()", excluded={0})
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
        return self._heat_capacity_gas_t(temperature)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._thermal_conductivity_liquid_tp(temperature, pressure)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._thermal_conductivity_gas_tp(temperature, pressure)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._viscosity_liquid_tp(temperature, pressure)

    @vectorize(signature="(),()->()", excluded={0})
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
        return self._viscosity_gas_tp(temperature, pressure)

    @vectorize(signature="(),()->()", excluded={0})
    def heat_capacity_solid_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:
        """Return the integral of solid heat capacity.

        Calculate the integral of solid heat capacity between temperature1
        and temperature2

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

    @vectorize(signature="(),()->()", excluded={0})
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

    @vectorize(signature="(),()->()", excluded={0})
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
