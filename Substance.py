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
    h_formation : float
        Standard state molar enthalpy of formation [J/mol], by default 
        None.
    h_formation_ig : float
        Ideal-gas molar enthalpy of formation [J/mol], by default None
    g_formation : float
        Standard state molar change of Gibbs energy of formation [J/mol]
        , by default None
    g_formation_ig : float
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
            self, name=None, mw=None, tc=None, pc=None, 
            omega=None, h_formation=None, h_formation_ig=None, 
            g_formation=None, g_formation_ig=None, volume_s_t=None, 
            volume_l_tp=None, volume_g_tp=None, heat_capacity_s_t=None, 
            heat_capacity_l_t=None, heat_capacity_g_t=None,
            thermal_conductivity_l_tp=None, 
            thermal_conductivity_g_tp=None, 
            viscosity_l_tp=None, viscosity_g_tp=None
        ):
        
        #Pure compound properties:
        self.name = name
        self.mw = mw
        self.tc = tc
        self.pc = pc
        self.omega = omega
        self.h_formation = h_formation
        self.h_formation_ig = h_formation_ig
        self.g_formation = g_formation
        self.g_formation_ig = g_formation_ig
        #Temperature dependent properties calculation functions:
        self._volume_s_t = volume_s_t
        self._volume_l_tp = volume_l_tp
        self._volume_g_tp = volume_g_tp
        self._heat_capacity_s_t = heat_capacity_s_t
        self._heat_capacity_l_t = heat_capacity_l_t
        self._heat_capacity_g_t = heat_capacity_g_t
        #self._thermal_conductivity_s_tp = thermal_conductivity_s_tp 
        self._thermal_conductivity_l_tp = thermal_conductivity_l_tp
        self._thermal_conductivity_g_tp = thermal_conductivity_g_tp
        self._viscosity_l_tp = viscosity_l_tp
        self._viscosity_g_tp = viscosity_g_tp
            
    @classmethod
    def from_thermo_database(cls, ID):
        """Method that use Bell Caleb's thermo library to construct the
        Substance object. Cite: Caleb Bell and Contributors (2016-2021).
        Thermo: Chemical properties component of Chemical Engineering 
        Design Library (ChEDL) https://github.com/CalebBell/thermo.

        Parameters
        ----------
        ID : string
            Name or CAS number of the substance
        """
        chemobj = Chemical(ID)

        substance_object = cls(
            name=chemobj.name, 
            mw=chemobj.MW, 
            tc=chemobj.Tc, 
            pc=chemobj.Pc, 
            omega=chemobj.omega, 
            h_formation=chemobj.Hfm, 
            h_formation_ig=chemobj.Hfgm, 
            g_formation=chemobj.Gfm, 
            g_formation_ig=chemobj.Gfgm, 
            volume_s_t=chemobj.VolumeSolid, 
            volume_l_tp=chemobj.VolumeLiquid, 
            volume_g_tp=chemobj.VolumeGas, 
            heat_capacity_s_t=chemobj.HeatCapacitySolid, 
            heat_capacity_l_t=chemobj.HeatCapacityLiquid, 
            heat_capacity_g_t=chemobj.HeatCapacityGas,
            thermal_conductivity_l_tp=chemobj.ThermalConductivityLiquid, 
            thermal_conductivity_g_tp=chemobj.ThermalConductivityGas, 
            viscosity_l_tp=chemobj.ViscosityLiquid, 
            viscosity_g_tp=chemobj.ViscosityGas
        )
        return substance_object

    def volume_solid(self, T):
        return self._volume_s_t(T)

    def volume_liquid(self, T, P):
        return self._volume_l_tp(T, P)

    def volume_gas(self, T, P):
        return self._volume_g_tp(T, P)

    def heat_capacity_solid(self, T):
        return self._heat_capacity_s_t(T)

    def heat_capacity_liquid(self, T):
        return self._heat_capacity_l_t(T)

    def heat_capacity_gas(self, T):
        return self._heat_capacity_g_t(T)

    #def thermal_conductivity_solid(self, T, P):
        #Not implemented (not in thermo library)
    #    return self._thermal_conductivity_s_tp(T, P)

    def thermal_conductivity_liquid(self, T, P):
        return self._thermal_conductivity_l_tp(T, P)  

    def thermal_conductivity_gas(self, T, P):
        return self._thermal_conductivity_g_tp(T, P)

    def viscosity_liquid(self, T, P):
        return self._viscosity_l_tp(T, P)

    def viscosity_gas(self, T, P):
        return self._viscosity_g_tp(T, P)