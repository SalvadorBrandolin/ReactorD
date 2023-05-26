import numpy as np

import pytest

import reactord as rd


def test_herning_zipperer():
    def pure_v1(t, p):
        return 3 * t + 2 * p

    def pure_v2(t, p):
        return 5 * t + 3 * p

    def pure_v3(t, p):
        return 10 * t + 4 * p

    a = rd.Substance(
        "A",
        molecular_weight=16,
        viscosity_liquid=pure_v1,
        viscosity_gas=pure_v1,
    )
    b = rd.Substance(
        "B",
        molecular_weight=9,
        viscosity_liquid=pure_v2,
        viscosity_gas=pure_v2,
    )
    c = rd.Substance(
        "C",
        molecular_weight=36,
        viscosity_liquid=pure_v3,
        viscosity_gas=pure_v3,
    )

    gas = rd.mix.IdealGas([a, b, c], viscosity_mixing_rule="grunberg_nissan")
    solution = rd.mix.IdealSolution(
        [a, b, c], viscosity_mixing_rule="grunberg_nissan"
    )

    temperatures = np.linspace(300, 400, 100)
    pressures = np.linspace(100000, 300000, 100)
    moles = np.array(
        [
            np.linspace(10, 80, 100),
            np.linspace(40, 60, 100),
            np.linspace(15, 150, 100),
        ]
    )

    mole_fractions = gas.mole_fractions(moles)

    assert (
        gas.mix_viscosity(mole_fractions, temperatures, pressures)
        == solution.mix_viscosity(mole_fractions, temperatures, pressures)
    ).all()

    pure_v = (
        (3 * temperatures + 2 * pressures),
        (5 * temperatures + 3 * pressures),
        (10 * temperatures + 4 * pressures),
    )

    mix_v_by_hand = np.exp(
        np.multiply(mole_fractions, np.log(pure_v)).sum(axis=0)
    )

    assert np.allclose(
        gas.mix_viscosity(mole_fractions, temperatures, pressures),
        solution.mix_viscosity(mole_fractions, temperatures, pressures),
    )

    assert np.allclose(
        gas.mix_viscosity(mole_fractions, temperatures, pressures),
        mix_v_by_hand,
    )

    gas.phase_nature = "A miserable little pile of secrets."

    with pytest.raises(ValueError):
        gas.mix_viscosity(mole_fractions, temperatures, pressures)
