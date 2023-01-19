import numpy as np

import pytest

import reactord as rd

from scipy.integrate import quad, quadrature

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
    method = chemical_obj.HeatCapacitySolids[0].method
    for t1, t2 in zip(temperature1, temperature2):
        for p in pressure:
            assert substance.heat_capacity_solid_dt_integral(
                t1, t2, p
            ) == chemical_obj.HeatCapacitySolids[0].calculate_integral(
                t1, t2, method
            )

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
    pressure = [101325, 200000, 300000, 400000]
    method = chemical_obj.HeatCapacityLiquids[0].method
    for t1, t2 in zip(temperature1, temperature2):
        for p in pressure:
            assert np.allclose(
                substance.heat_capacity_liquid_dt_integral(t1, t2, p),
                (
                    chemical_obj.HeatCapacityLiquids[0].calculate_integral(
                        t1, t2, method
                    )
                ),
            )

            assert np.allclose(
                substance.heat_capacity_liquid_dt_integral(t1, t2, p),
                quad(
                    chemical_obj.HeatCapacityLiquids[0].T_dependent_property,
                    t1,
                    t2,
                )[0],
                atol=1000
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
                (
                    chemical_obj.HeatCapacityGases[0].calculate_integral(
                        t1, t2, method
                    )
                ),
            )

            assert np.allclose(
                substance.heat_capacity_gas_dt_integral(t1, t2, p),
                quad(
                    chemical_obj.HeatCapacityGases[0].T_dependent_property,
                    t1,
                    t2,
                )[0],
            )
