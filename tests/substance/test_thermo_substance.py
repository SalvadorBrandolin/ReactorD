import numpy as np

import pytest

import reactord as rd

from scipy.integrate import quad

from thermo import ChemicalConstantsPackage


compounds = [("water"), ("ethanol"), ("pentane"), ("toluene")]


@pytest.mark.parametrize("name", compounds)
def test_name(name):
    substance = rd.Substance.from_thermo_database("some_name", name)
    assert substance.name == "some_name"


@pytest.mark.parametrize("name", compounds)
def test_molecular_weight(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.molecular_weight == chemical_obj.constants.MWs[0]


@pytest.mark.parametrize("name", compounds)
def test_normal_boiling_point(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.normal_boiling_point == chemical_obj.constants.Tbs[0]


@pytest.mark.parametrize("name", compounds)
def test_normal_melting_point(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.normal_melting_point == chemical_obj.constants.Tms[0]


@pytest.mark.parametrize("name", compounds)
def test_critical_temperature(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.critical_temperature == chemical_obj.constants.Tcs[0]


@pytest.mark.parametrize("name", compounds)
def test_critical_pressure(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.critical_pressure == chemical_obj.constants.Pcs[0]


@pytest.mark.parametrize("name", compounds)
def test_acentric_factor(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.acentric_factor == chemical_obj.constants.omegas[0]


@pytest.mark.parametrize("name", compounds)
def test_formation_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.formation_enthalpy == chemical_obj.constants.Hf_STPs[0]


@pytest.mark.parametrize("name", compounds)
def test_formation_enthalpy_ig(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.formation_enthalpy_ig == chemical_obj.constants.Hfgs[0]


@pytest.mark.parametrize("name", compounds)
def test_formation_gibbs_ig(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    assert substance.formation_gibbs_ig == chemical_obj.constants.Gfgs[0]


@pytest.mark.parametrize("name", compounds)
def test_vectorize_functions(name):
    substance = rd.Substance.from_thermo_database(name, name)
    assert substance.vectorize_functions == True


@pytest.mark.parametrize("name", compounds)
def test_vaporization_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.vaporization_enthalpy(t) == (
            chemical_obj.EnthalpyVaporizations[0].T_dependent_property(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_sublimation_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.sublimation_enthalpy(t) == (
            chemical_obj.EnthalpySublimations[0].T_dependent_property(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_fusion_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 350, 400]
    for t in temperature:
        assert substance.fusion_enthalpy(t) == (
            chemical_obj.EnthalpySublimations[0].T_dependent_property(t)
            - chemical_obj.EnthalpyVaporizations[0].T_dependent_property(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_volume_solid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.volume_solid(t, p) == (
                chemical_obj.VolumeSolids[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_volume_liquid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.volume_liquid(t, p) == (
                chemical_obj.VolumeLiquids[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_volume_gas(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            method = chemical_obj.VolumeGases[0].method_P
            assert substance.volume_gas(t, p) == (
                chemical_obj.VolumeGases[0].calculate_P(t, p, method)
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_solid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.heat_capacity_solid(t, p) == (
                chemical_obj.HeatCapacitySolids[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_liquid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.heat_capacity_liquid(t, p) == (
                chemical_obj.HeatCapacityLiquids[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_gas(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.heat_capacity_gas(t, 101325) == (
                chemical_obj.HeatCapacityGases[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_thermal_conductivity_liquid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.thermal_conductivity_liquid(t, p) == (
                chemical_obj.ThermalConductivityLiquids[
                    0
                ].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_thermal_conductivity_gas(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.thermal_conductivity_gas(t, p) == (
                chemical_obj.ThermalConductivityGases[0].T_dependent_property(
                    t
                )
            )


@pytest.mark.parametrize("name", compounds)
def test_viscosity_liquid(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.viscosity_liquid(t, p) == (
                chemical_obj.ViscosityLiquids[0].T_dependent_property(t)
            )


@pytest.mark.parametrize("name", compounds)
def test_viscosity_gas(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert np.allclose(
                substance.viscosity_gas(t, p),
                (chemical_obj.ViscosityGases[0].T_dependent_property(t)),
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_solid_dt_integral(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature1 = [300, 400, 500, 600, 700]
    temperature2 = [400, 500, 600, 700, 800]
    pressure = [101325, 200000, 300000, 400000]
    for t1, t2 in zip(temperature1, temperature2):
        for p in pressure:
            assert np.allclose(
                substance.heat_capacity_solid_dt_integral(t1, t2, p),
                quad(
                    chemical_obj.HeatCapacitySolids[0].T_dependent_property,
                    t1,
                    t2,
                )[0],
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_liquid_dt_integral(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature1 = [300, 400, 500, 600, 700]
    temperature2 = [400, 500, 600, 700, 800]
    pressure = [101325, 202325, 303325]
    method = chemical_obj.HeatCapacityLiquids[0].method
    for t1, t2 in zip(temperature1, temperature2):
        for p in pressure:
            assert np.allclose(
                substance.heat_capacity_liquid_dt_integral(t1, t2, p),
                quad(
                    chemical_obj.HeatCapacityLiquids[0].T_dependent_property,
                    t1,
                    t2,
                )[0],
            )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_gas_dt_integral(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature1 = [300, 400, 500, 600, 700]
    temperature2 = [400, 500, 600, 700, 800]
    pressure = [101325, 200000, 300000, 400000]
    method = chemical_obj.HeatCapacityGases[0].method
    for t1, t2 in zip(temperature1, temperature2):
        for p in pressure:
            assert np.allclose(
                substance.heat_capacity_gas_dt_integral(t1, t2, p),
                quad(
                    chemical_obj.HeatCapacityGases[0].T_dependent_property,
                    t1,
                    t2,
                )[0],
            )

@pytest.mark.parametrize("name", compounds)
def test_all_together_vectorized(name):
    substance = rd.Substance.from_thermo_database(name, name)
    chemical_obj = ChemicalConstantsPackage.correlations_from_IDs([name])
    temperature1 = np.array([300, 350, 400, 450])
    temperature2 = np.array([300, 350, 400, 450])
    pressure = np.array([101325, 201325, 301325, 404325])

    # Constants
    assert substance.name == name
    assert substance.molecular_weight == chemical_obj.constants.MWs[0]
    assert substance.normal_boiling_point == chemical_obj.constants.Tbs[0]
    assert substance.normal_melting_point == chemical_obj.constants.Tms[0]
    assert substance.critical_temperature == chemical_obj.constants.Tcs[0]
    assert substance.critical_pressure == chemical_obj.constants.Pcs[0]
    assert substance.acentric_factor == chemical_obj.constants.omegas[0]
    assert substance.formation_enthalpy == chemical_obj.constants.Hf_STPs[0]
    assert substance.formation_enthalpy_ig == chemical_obj.constants.Hfgs[0]
    assert substance.formation_gibbs_ig == chemical_obj.constants.Gfgs[0]
    assert substance.vectorize_functions == True

    # Functions
    assert isinstance(substance._vaporization_enthalpy, np.vectorize)
    assert isinstance(substance._sublimation_enthalpy, np.vectorize)
    assert isinstance(substance._volume_solid, np.vectorize)
    assert isinstance(substance._volume_liquid, np.vectorize)
    assert isinstance(substance._volume_gas, np.vectorize)
    assert isinstance(substance._heat_capacity_solid, np.vectorize)
    assert isinstance(substance._heat_capacity_liquid, np.vectorize)
    assert isinstance(substance._heat_capacity_gas, np.vectorize)
    assert isinstance(substance._thermal_conductivity_liquid, np.vectorize)
    assert isinstance(substance._thermal_conductivity_gas, np.vectorize)
    assert isinstance(substance._viscosity_liquid, np.vectorize)
    assert isinstance(substance._viscosity_gas, np.vectorize)
    assert isinstance(substance._heat_capacity_solid_dt_integral, np.vectorize)
    assert isinstance(
        substance._heat_capacity_liquid_dt_integral, np.vectorize
    )
    assert isinstance(substance._heat_capacity_gas_dt_integral, np.vectorize)

    #Eval
    reactord_vaporization_enthalpy = substance.vaporization_enthalpy(temperature1)
    reactord_sublimation_enthalpy = substance.sublimation_enthalpy(temperature1)
    reactord_fusion_enthalpy = substance.fusion_enthalpy(temperature1)
    reactord_volume_solid = substance.volume_solid(temperature1, pressure)
    reactord_volume_liquid = substance.volume_liquid(temperature1, pressure)
    reactord_volume_gas = substance.volume_gas(temperature1, pressure)
    reactord_heat_capacity_solid = substance.heat_capacity_solid(temperature1, pressure)
    reactord_heat_capacity_liquid = substance.heat_capacity_liquid(temperature1, pressure)
    reactord_heat_capacity_gas = substance.heat_capacity_gas(temperature1, pressure)
    reactord_thermal_conductivity_liquid = substance.thermal_conductivity_liquid(temperature1, pressure)
    reactord_thermal_conductivity_gas = substance.thermal_conductivity_gas(temperature1, pressure)
    reactord_viscosity_liquid = substance.viscosity_liquid(temperature1, pressure)
    reactord_viscosity_gas = substance.viscosity_gas(temperature1, pressure)
    reactord_heat_capacity_solid_dt_integral = substance.heat_capacity_solid_dt_integral(temperature1, temperature2, pressure)
    reactord_heat_capacity_liquid_dt_integral = substance.heat_capacity_liquid_dt_integral(temperature1, temperature2, pressure)
    reactord_heat_capacity_gas_dt_integral = substance.heat_capacity_gas_dt_integral(temperature1, temperature2, pressure)
