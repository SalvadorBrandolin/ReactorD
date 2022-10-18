from thermo.chemical import Chemical

import reactord as rd


def test_equal_evaluation():
    substance = rd.Substance.from_thermo_database("methane")
    chemical_obj = Chemical("methane")

    assert substance.name == chemical_obj.name
    assert substance.mw == chemical_obj.MW
    assert substance.normal_boiling_point == chemical_obj.Tb
    assert substance.normal_melting_point == chemical_obj.Tm
    assert substance.tc == chemical_obj.Tc
    assert substance.pc == chemical_obj.Pc
    assert substance.omega == chemical_obj.omega
    assert substance.formation_enthalpy == chemical_obj.Hfm
    assert substance.formation_enthalpy_ig == chemical_obj.Hfgm

    temperature = [300, 400, 500, 600, 700]
    pressure = [101325, 200000, 300000, 400000]

    for t in temperature:
        assert substance.vaporization_enthalpy(t) == (
            chemical_obj.EnthalpyVaporization(t)
        )

        assert substance.sublimation_enthalpy(t) == (
            chemical_obj.EnthalpySublimation(t)
        )

        assert substance.volume_solid(t) == (chemical_obj.VolumeSolid(t))

        assert substance.heat_capacity_solid(t) == (
            chemical_obj.HeatCapacitySolid(t)
        )

        assert substance.heat_capacity_liquid(t) == (
            chemical_obj.HeatCapacityLiquid(t)
        )

        assert substance.heat_capacity_gas(t) == (
            chemical_obj.HeatCapacityGas(t)
        )

        for p in pressure:
            assert substance.volume_liquid(t, p) == (
                chemical_obj.VolumeLiquid(t, p)
            )

            assert substance.volume_gas(t, p) == (chemical_obj.VolumeGas(t, p))

            assert substance.thermal_conductivity_liquid(t, p) == (
                chemical_obj.ThermalConductivityLiquid(t, p)
            )

            assert substance.thermal_conductivity_gas(t, p) == (
                chemical_obj.ThermalConductivityGas(t, p)
            )

            assert substance.viscosity_liquid(t, p) == (
                chemical_obj.ViscosityLiquid(t, p)
            )

            assert substance.viscosity_gas(t, p) == (
                chemical_obj.ViscosityGas(t, p)
            )
