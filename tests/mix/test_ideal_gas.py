import numpy as np

import reactord as rd


def test_one_substance_mix():
    methane = rd.Substance.from_thermo_database("methane")

    oxygen = rd.Substance.from_thermo_database("oxygen")

    hydrogen = rd.Substance.from_thermo_database("hydrogen")

    mix1 = rd.mix.IdealGas(A=methane)
    mix2 = rd.mix.IdealGas(B=oxygen)
    mix3 = rd.mix.IdealGas(C=hydrogen)

    r = 8.31446261815324  # m3⋅Pa/K/mol

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000])

    # ------------------- Tests Adrian ----------------------------------------
    # Test of __str__ method: --Tests Adrian--
    assert (
        mix1.__str__()
        == "The mixture contains the following 1 components:\nMethane\n"
    )
    assert (
        mix2.__str__()
        == "The mixture contains the following 1 components:\nOxygen\n"
    )
    assert (
        mix3.__str__()
        == "The mixture contains the following 1 components:\nHydrogen\n"
    )

    assert mix1._formation_enthalpies_set() == methane.formation_enthalpy_ig
    assert mix2._formation_enthalpies_set() == oxygen.formation_enthalpy_ig
    assert mix3._formation_enthalpies_set() == hydrogen.formation_enthalpy_ig

    for t in temperature:
        assert mix1.formation_enthalpies_correction(
            t
        ) == methane.heat_capacity_gas_dt_integral(298.15, t)
        assert mix2.formation_enthalpies_correction(
            t
        ) == oxygen.heat_capacity_gas_dt_integral(298.15, t)
        assert mix3.formation_enthalpies_correction(
            t
        ) == hydrogen.heat_capacity_gas_dt_integral(298.15, t)
    # ----------------------------------------------------------------------

    for n in compositions:
        assert mix1.mol_fractions(n) == 1.0
        assert mix2.mol_fractions(n) == 1.0
        assert mix3.mol_fractions(n) == 1.0

        for t, p in zip(temperature, pressure):
            assert mix1.concentrations(n, t, p) == p / (r * t)
            assert mix1.concentrations(n, t, p) == mix2.concentrations(n, t, p)
            assert mix2.concentrations(n, t, p) == mix3.concentrations(n, t, p)

            assert n * r * t / p == mix2.volume(n, t, p)
            assert mix2.volume(n, t, p) == mix3.volume(n, t, p)
            assert mix1.volume(n, t, p) == mix2.volume(n, t, p)

            assert mix1.mix_heat_capacity(n, t, p) == (
                methane.heat_capacity_gas(t)
            )

            assert mix2.mix_heat_capacity(n, t, p) == (
                oxygen.heat_capacity_gas(t)
            )

            assert mix3.mix_heat_capacity(n, t, p) == (
                hydrogen.heat_capacity_gas(t)
            )

            # Partial pressures of AbstractMix tested: --Tests Adrian--
            assert mix1.partial_pressures(n, t, p) == p
            assert mix2.partial_pressures(n, t, p) == p
            assert mix3.partial_pressures(n, t, p) == p


def test_three_substances_mix():
    co2 = rd.Substance.from_thermo_database("co2")
    ethane = rd.Substance.from_thermo_database("ethane")
    chlorine = rd.Substance.from_thermo_database("chlorine")

    mixture = rd.mix.IdealGas(A=co2, B=ethane, C=chlorine)

    r = 8.31446261815324  # m3⋅Pa/K/mol

    compositions = np.array(
        [
            [1, 5, 8],
            [10, 15, 20],
            [100, 50, 30],
            [1000, 8000, 500],
            [10000, 500, 4000],
        ]
    )

    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000])

    for t, p in zip(temperature, pressure):
        raw_formation_enthalpies_correction = np.array(
            [
                co2.heat_capacity_gas_dt_integral(298.15, t),
                ethane.heat_capacity_gas_dt_integral(298.15, t),
                chlorine.heat_capacity_gas_dt_integral(298.15, t),
            ]
        )

        raw_heat_capacities = np.array(
            [
                co2.heat_capacity_gas(t),
                ethane.heat_capacity_gas(t),
                chlorine.heat_capacity_gas(t),
            ]
        )

        raw_densities = p / (r * t)

        for moles in compositions:
            raw_mol_fractions = np.divide(moles, np.sum(moles, axis=0))
            raw_concentrations = np.multiply(raw_mol_fractions, raw_densities)

            # Test of concentrations method
            assert (
                raw_concentrations == mixture.concentrations(moles, t, p)
            ).all()  # OKAY

            # Test of volume method
            raw_vol = sum(moles) * r * t / p
            assert mixture.volume(moles, t, p) == raw_vol  # OKAY

            # Test of mix_heat_capacity method
            raw_mix_heat_capacity = np.dot(
                raw_heat_capacities, raw_mol_fractions
            )
            assert raw_mix_heat_capacity == mixture.mix_heat_capacity(
                moles, t, p
            )  # OKAY

            # Test of formation_enthalpies_set method
            raw_enthalpies_set = np.array(
                [
                    co2.formation_enthalpy_ig,
                    ethane.formation_enthalpy_ig,
                    chlorine.formation_enthalpy_ig,
                ]
            )

            assert (
                mixture._formation_enthalpies_set() == raw_enthalpies_set
            ).all()  # OKAY

            # Test of formation_enthalpies_correction method
            assert (
                raw_formation_enthalpies_correction
                == mixture.formation_enthalpies_correction(t)
            ).all()  # OKAY
