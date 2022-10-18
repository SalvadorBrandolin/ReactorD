from scipy.integrate import quad
from thermo.chemical import Chemical


class Substance:
    """Substance object class

    Parameters
    ----------
    name : string
        Name of the substance, by default None
    mw : float
        Molecular weigth of the substance [g/mol], by default None
    tc : float
        Critical temperature of the substance [K], by default None
    pc : float
        Critical pressure of the substance [Pa], by default None
    omega : float
        Acentric factor of the substance, by default None
    formation_enthalpy : float
        Standard state molar enthalpy of formation [J/mol], by default
        None.
    formation_enthalpy_ig : float
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    formation_gibbs : float
        Standard state molar change of Gibbs energy of formation [J/mol]
        , by default None
    formation_gibbs_ig : float
        Ideal-gas molar change of Gibbs energy of formation [J/mol], by
        default None
    volume_s_t : function
        Solid molar volume as a python function of temperature
        volume_solid(T) [m^3/mol], by default None
    volume_l_tp : function
        Liquid molar volume as a python function of temperature [K]
        and pressure [Pa], volume_liquid(T, P) [m^3/mol], by default
        None
    volume_g_tp : function
        Gas molar volume as a python function of temperature [K]
        and pressure [Pa], volume_gas(T, P) [m^3/mol], by default None
    heat_capacity_s_t : function
        Solid heat capacity as a python function of temperature
        [K], heat_capacity_solid(T) [J/mol/K], by default None
    heat_capacity_l_t : function
        Liquid heat capacity as a python function of temperature
        [K], heat_capacity_liquid(T) [J/mol/K], by default None
    heat_capacity_g_t : function
        Ideal-gas heat capacity as a python function of temperature
        [K], heat_capacity_gas(T) [J/mol/K], by default None
    thermal_conductivity_l_tp : function
        Liquid thermal conductivity as a python function of temperature
        [K] and pressure [Pa], thermal_conductivity_liquid(T) [W/m/K],
        by default None
    thermal_conductivity_g_tp : function
        Gas thermal conductivity as a python function of temperature [K]
        and pressure [Pa], thermal_conductivity_liquid(T) [W/m/K], by
        default None
    viscosity_l_tp : function
        Liquid viscocity as a python function of temperature [K] and
        pressure [Pa], viscosity_liquid(T) [W/m/K], by default None
    viscosity_g_tp : function
        Gas viscocity as a python function of temperature [K] and
        pressure [Pa], viscosity_gas(T) [W/m/K], by default None
    """

    def __init__(
        self,
        name=None,
        mw=None,
        normal_boiling_point=None,
        normal_melting_point=None,
        tc=None,
        pc=None,
        omega=None,
        formation_enthalpy=None,
        formation_enthalpy_ig=None,
        formation_gibbs=None,
        formation_gibbs_ig=None,
        vaporization_enthalpy_t=None,
        sublimation_enthalpy_t=None,
        volume_s_t=None,
        volume_l_tp=None,
        volume_g_tp=None,
        heat_capacity_s_t=None,
        heat_capacity_l_t=None,
        heat_capacity_g_t=None,
        thermal_conductivity_l_tp=None,
        thermal_conductivity_g_tp=None,
        viscosity_l_tp=None,
        viscosity_g_tp=None,
    ):

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
        # Temperature dependent properties calculation functions:
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
    def from_thermo_database(cls, identificator):
        """Method that use Bell Caleb's thermo library to construct the
        Substance object. Cite: Caleb Bell and Contributors (2016-2021).
        Thermo: Chemical properties component of Chemical Engineering
        Design Library (ChEDL) https://github.com/CalebBell/thermo.

        Parameters
        ----------
        identificator : string
            Name or CAS number of the substance
        """
        chemobj = Chemical(identificator)

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

    def vaporization_enthalpy(self, temperature):
        return self._vaporization_enthalpy_t(temperature)

    def sublimation_enthalpy(self, temperature):
        return self._sublimation_enthalpy_t(temperature)

    def fusion_enthalpy(self, temperature):
        fusion_h = self._sublimation_enthalpy_t(
            temperature
        ) - self._vaporization_enthalpy_t(temperature)
        return fusion_h

    def volume_solid(self, temperature):
        return self._volume_s_t(temperature)

    def volume_liquid(self, temperature, pressure):
        return self._volume_l_tp(temperature, pressure)

    def volume_gas(self, temperature, pressure):
        return self._volume_g_tp(temperature, pressure)

    def heat_capacity_solid(self, temperature):
        return self._heat_capacity_s_t(temperature)

    def heat_capacity_liquid(self, temperature):
        return self._heat_capacity_l_t(temperature)

    def heat_capacity_gas(self, temperature):
        return self._heat_capacity_g_t(temperature)

    # def thermal_conductivity_solid(self, temperature, pressure):
    # Not implemented (not in thermo library)
    #    return self._thermal_conductivity_s_tp(temperature, pressure)

    def thermal_conductivity_liquid(self, temperature, pressure):
        return self._thermal_conductivity_l_tp(temperature, pressure)

    def thermal_conductivity_gas(self, temperature, pressure):
        return self._thermal_conductivity_g_tp(temperature, pressure)

    def viscosity_liquid(self, temperature, pressure):
        return self._viscosity_l_tp(temperature, pressure)

    def viscosity_gas(self, temperature, pressure):
        return self._viscosity_g_tp(temperature, pressure)

    def heat_capacity_solid_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:

        integral, err = quad(
            self.heat_capacity_solid, a=temperature1, b=temperature2
        )

        return integral

    def heat_capacity_liquid_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:

        integral, err = quad(
            self.heat_capacity_liquid, a=temperature1, b=temperature2
        )

        return integral

    def heat_capacity_gas_dt_integral(
        self, temperature1: float, temperature2: float
    ) -> float:

        integral, err = quad(
            self.heat_capacity_gas, a=temperature1, b=temperature2
        )

        return integral
