"""Thermo substance constructor module.

Function to define a substance from Bell Calleb Thermo database.Cite:
Caleb Bell and Contributors (2016-2021). Thermo: Chemical properties component
of Chemical Engineering Design Library (ChEDL)
https://github.com/CalebBell/thermo.
"""

from scipy.integrate import quad

from thermo import ChemicalConstantsPackage


def thermo_substance_constructor(cls, name: str, thermo_identification: str):
    """Substance constructor from Thermo database.

    The functions of the Thermo library are all exclusively temperature
    functions (pressure is not used), except the volume_gas method which it
    uses both temperature and pressure.

    Parameters
    ----------
    cls : Substance
        Substance class.
    name : str
        Name that will be assigned to the Substance object.
    thermo_identification : str
        Name or CAS number of the substance that will be used to search in the
        Thermo library.

    Returns
    -------
    Substance
        Instantiated Substance object from thermo database.
    """
    corr = ChemicalConstantsPackage.correlations_from_IDs(
        [thermo_identification]
    )

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

    def heat_capacity_solid(temperature: float, pressure: float) -> float:
        heat_cap = corr.HeatCapacitySolids[0].T_dependent_property(temperature)
        return heat_cap

    def heat_capacity_liquid(temperature: float, pressure: float) -> float:
        heat_cap = corr.HeatCapacityLiquids[0].T_dependent_property(
            temperature
        )
        return heat_cap

    def heat_capacity_gas(temperature: float, pressure: float) -> float:
        heat_cap = corr.HeatCapacityGases[0].T_dependent_property(temperature)
        return heat_cap

    def thermal_conductivity_liquid(
        temperature: float, pressure: float
    ) -> float:
        thermal_cond = corr.ThermalConductivityLiquids[0].T_dependent_property(
            temperature
        )
        return thermal_cond

    def thermal_conductivity_gas(temperature: float, pressure: float) -> float:
        thermal_cond = corr.ThermalConductivityGases[0].T_dependent_property(
            temperature
        )
        return thermal_cond

    def viscosity_liquid(temperature: float, pressure: float) -> float:
        viscosity_liquid = corr.ViscosityLiquids[0].T_dependent_property(
            temperature
        )
        return viscosity_liquid

    def viscosity_gas(temperature: float, pressure: float) -> float:
        viscosity_gas = corr.ViscosityGases[0].T_dependent_property(
            temperature
        )
        return viscosity_gas

    def heat_capacity_solid_dt_integral(
        temperature1: float, temperature2: float, pressure: float
    ) -> float:
        integral = quad(
            corr.HeatCapacitySolids[0].T_dependent_property,
            temperature1,
            temperature2,
        )[0]
        return integral

    def heat_capacity_liquid_dt_integral(
        temperature1: float, temperature2: float, pressure: float
    ) -> float:
        integral = quad(
            corr.HeatCapacityLiquids[0].T_dependent_property,
            temperature1,
            temperature2,
        )[0]
        return integral

    def heat_capacity_gas_dt_integral(
        temperature1: float, temperature2: float, pressure: float
    ) -> float:
        integral = quad(
            corr.HeatCapacityGases[0].T_dependent_property,
            temperature1,
            temperature2,
        )[0]
        return integral

    substance_object = cls(
        name=name,
        molecular_weight=corr.constants.MWs[0],
        normal_boiling_point=corr.constants.Tbs[0],
        normal_melting_point=corr.constants.Tms[0],
        critical_temperature=corr.constants.Tcs[0],
        critical_pressure=corr.constants.Pcs[0],
        acentric_factor=corr.constants.omegas[0],
        formation_enthalpy=corr.constants.Hf_STPs[0],
        formation_enthalpy_ig=corr.constants.Hfgs[0],
        vaporization_enthalpy=vaporization_enthalpy,
        sublimation_enthalpy=sublimation_enthalpy,
        volume_solid=volume_solid,
        volume_liquid=volume_liquid,
        heat_capacity_solid=heat_capacity_solid,
        heat_capacity_liquid=heat_capacity_liquid,
        heat_capacity_gas=heat_capacity_gas,
        thermal_conductivity_liquid=thermal_conductivity_liquid,
        thermal_conductivity_gas=thermal_conductivity_gas,
        viscosity_liquid=viscosity_liquid,
        viscosity_gas=viscosity_gas,
        heat_capacity_solid_dt_integral=heat_capacity_solid_dt_integral,
        heat_capacity_liquid_dt_integral=heat_capacity_liquid_dt_integral,
        heat_capacity_gas_dt_integral=heat_capacity_gas_dt_integral,
        vectorize_functions=True,
    )

    return substance_object
