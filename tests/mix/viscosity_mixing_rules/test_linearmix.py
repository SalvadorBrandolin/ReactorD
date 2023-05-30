import numpy as np

import pytest

import reactord as rd


def test_linear():
    def pure_v1(t, p):
        return 3 * t + 2 * p

    def pure_v2(t, p):
        return 5 * t + 3 * p

    def pure_v3(t, p):
        return 10 * t + 4 * p

    a = rd.Substance("A", viscosity_liquid=pure_v1, viscosity_gas=pure_v1)
    b = rd.Substance("B", viscosity_liquid=pure_v2, viscosity_gas=pure_v2)
    c = rd.Substance("C", viscosity_liquid=pure_v3, viscosity_gas=pure_v3)

    gas = rd.mix.IdealGas([a, b, c], viscosity_mixing_rule="linear")
    solution = rd.mix.IdealSolution([a, b, c], viscosity_mixing_rule="linear")

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

    mix_v_by_hand = (
        (mole_fractions[0] * (3 * temperatures + 2 * pressures))
        + (mole_fractions[1] * (5 * temperatures + 3 * pressures))
        + (mole_fractions[2] * (10 * temperatures + 4 * pressures))
    )

    assert (
        gas.mix_viscosity(mole_fractions, temperatures, pressures)
        == solution.mix_viscosity(mole_fractions, temperatures, pressures)
    ).all()
    assert (
        gas.mix_viscosity(mole_fractions, temperatures, pressures)
        == mix_v_by_hand
    ).all()

    gas.phase_nature = "Kratos"

    with pytest.raises(ValueError):
        gas.mix_viscosity(mole_fractions, temperatures, pressures)
