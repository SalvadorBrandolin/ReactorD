import numpy as np

import reactord as rd


def test_one_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane")

    mixture = rd.mix.IdealSolution([hexane])

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = 101325

    for z in compositions:
        assert mixture.mol_fracations(z) == 1.0

        for t in temperature:
            assert mixture.volume(z, t, pressure) == (
                hexane.volume_liquid(t, pressure),
            )

            assert mixture.concentrations(z, t, pressure) == (
                (1 / hexane.volume_liquid(t, pressure))
            )

            assert mixture.mix_heat_capacity(z, t, pressure) == (
                hexane.heat_capacity_liquid(t)
            )
