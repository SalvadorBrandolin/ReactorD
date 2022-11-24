import pytest

import reactord as rd

from thermo.chemical import Chemical

compounds = [("water"), ("methane"), ("pentane"), ("toluene")]


@pytest.mark.parametrize("name", compounds)
def test_name(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.name == chemical_obj.name


@pytest.mark.parametrize("name", compounds)
def test_mw(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.mw == chemical_obj.MW


@pytest.mark.parametrize("name", compounds)
def test_normal_boiling_point(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.normal_boiling_point == chemical_obj.Tb


@pytest.mark.parametrize("name", compounds)
def test_normal_melting_point(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.normal_melting_point == chemical_obj.Tm


@pytest.mark.parametrize("name", compounds)
def test_tc(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.tc == chemical_obj.Tc


@pytest.mark.parametrize("name", compounds)
def test_pc(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.pc == chemical_obj.Pc


@pytest.mark.parametrize("name", compounds)
def test_omega(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.omega == chemical_obj.omega


@pytest.mark.parametrize("name", compounds)
def test_formation_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.formation_enthalpy == chemical_obj.Hfm


@pytest.mark.parametrize("name", compounds)
def test_formation_enthalpy_ig(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    assert substance.formation_enthalpy_ig == chemical_obj.Hfgm


@pytest.mark.parametrize("name", compounds)
def test_vaporization_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.vaporization_enthalpy(t) == (
            chemical_obj.EnthalpyVaporization(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_sublimation_enthalpy(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.sublimation_enthalpy(t) == (
            chemical_obj.EnthalpySublimation(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_volume_solid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.volume_solid(t) == (chemical_obj.VolumeSolid(t))


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_solid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.heat_capacity_solid(t) == (
            chemical_obj.HeatCapacitySolid(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_liquid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.heat_capacity_liquid(t) == (
            chemical_obj.HeatCapacityLiquid(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_heat_capacity_gas(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    for t in temperature:
        assert substance.heat_capacity_gas(t) == (
            chemical_obj.HeatCapacityGas(t)
        )


@pytest.mark.parametrize("name", compounds)
def test_volume_liquid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.volume_liquid(t, p) == (
                chemical_obj.VolumeLiquid(t, p)
            )


@pytest.mark.parametrize("name", compounds)
def test_volume_gas(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.volume_gas(t, p) == (chemical_obj.VolumeGas(t, p))


@pytest.mark.parametrize("name", compounds)
def test_thermal_conductivity_liquid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.thermal_conductivity_liquid(t, p) == (
                chemical_obj.ThermalConductivityLiquid(t, p)
            )


@pytest.mark.parametrize("name", compounds)
def test_thermal_conductivity_gas(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.thermal_conductivity_gas(t, p) == (
                chemical_obj.ThermalConductivityGas(t, p)
            )


@pytest.mark.parametrize("name", compounds)
def test_viscosity_liquid(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.viscosity_liquid(t, p) == (
                chemical_obj.ViscosityLiquid(t, p)
            )


@pytest.mark.parametrize("name", compounds)
def test_viscosity_gas(name):
    substance = rd.Substance.from_thermo_database(name)
    chemical_obj = Chemical(name)
    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]
    for t in temperature:
        for p in pressure:
            assert substance.viscosity_gas(t, p) == (
                chemical_obj.ViscosityGas(t, p)
            )
