import numpy as np

import reactord as rd


def test_one_substance_mix():
    methane = rd.Substance.from_thermo_database("methane")

    oxygen = rd.Substance.from_thermo_database("oxygen")

    hydrogen = rd.Substance.from_thermo_database("hydrogen")

    mix1 = rd.mix.IdealGas([methane])
    mix2 = rd.mix.IdealGas([oxygen])
    mix3 = rd.mix.IdealGas([hydrogen])

    r = 8.31446261815324  # m3⋅Pa/K/mol

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000])

    assert mix1._formation_enthalpies_set() == methane.formation_enthalpy_ig
    assert mix2._formation_enthalpies_set() == oxygen.formation_enthalpy_ig
    assert mix3._formation_enthalpies_set() == hydrogen.formation_enthalpy_ig
    
    #for t in temperature:
    #assert mix1.formation_enthalpies_correction(t) == methane.formation_entha
    

    for n in compositions:
        assert mix1.mol_fracations(n) == 1.0
        assert mix2.mol_fracations(n) == 1.0
        assert mix3.mol_fracations(n) == 1.0

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
