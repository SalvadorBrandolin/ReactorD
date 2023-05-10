import numpy as np

import reactord as rd
from reactord.kinetic.reaction_enthalpy import dh_not_specified, dh_specified

from scipy.integrate import quad

from thermo import ChemicalConstantsPackage


def test_dh_specified():
    a = rd.Substance("a")
    b = rd.Substance("b")
    c = rd.Substance("c")
    d = rd.Substance("d")

    mix = rd.mix.IdealGas([a, b, c, d])

    def rate():
        ...

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a + b > c + d, "rate": rate, "DH": -10_000},
            "r2": {"eq": a + 2 * b > 3 * c + d, "rate": rate, "DH": 20_000},
            "r3": {"eq": 2 * a + b > c + d, "rate": rate, "DH": -30_000},
        },
        kinetic_constants={},
    )

    temperature = np.linspace(200, 600, 20)
    pressure = np.linspace(101325, 301325, 20)

    kinetic.set_dh_function()

    # Compare one by one
    for t, p in zip(temperature, pressure):
        assert (
            np.array([-10000, 20000, -30000]) == dh_specified(kinetic, t, p).T
        ).all()

        assert (
            kinetic.dhs_evaluate(t, p) == dh_specified(kinetic, t, p)
        ).all()

    # Vectorial evaluation
    dhs = np.ones((3, 20))
    dhs[0, :] = np.full(20, -10_000)
    dhs[1, :] = np.full(20, 20_000)
    dhs[2, :] = np.full(20, -30_000)

    assert (kinetic.dhs_evaluate(temperature, pressure) == dhs).all()


def test_dh_not_specified_one_reaction():
    a = rd.Substance.from_thermo_database("ch4", "methane")
    b = rd.Substance.from_thermo_database("o2", "oxygen")
    c = rd.Substance.from_thermo_database("co2", "carbon dioxide")
    d = rd.Substance.from_thermo_database("h2o", "water")

    mix = rd.mix.IdealGas([a, b, c, d])

    def rate():
        ...

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"combustion": {"eq": a + 2 * b > c + 2 * d, "rate": rate}},
        kinetic_constants={},
    )

    kinetic.set_dh_function()

    cons, corr = ChemicalConstantsPackage.from_IDs(
        ["methane", "oxygen", "carbon dioxide", "water"]
    )

    hs = cons.Hfgs
    cps = corr.HeatCapacityGases

    deltah_std = (hs[2] + 2 * hs[3]) - (hs[0] + 2 * hs[1])

    assert kinetic.dhs_evaluate(298.15, 101325) == deltah_std

    temperatures = np.linspace(298.15, 600, 100)
    pressures = np.full(100, 101325)

    dhs_store = np.array([])
    # One by one
    for t in temperatures:
        integrals = np.array(
            [quad(cp.T_dependent_property, 298.15, t)[0] for cp in cps]
        )

        dh_corrected = (
            deltah_std
            + (integrals[2] + 2 * integrals[3])
            - (integrals[0] + 2 * integrals[1])
        )

        assert np.allclose(
            dh_corrected, kinetic.dhs_evaluate(t, 101325), atol=1e-10
        )

        dhs_store = np.append(dhs_store, dh_corrected)

    # vectorial call
    assert np.allclose(
        dhs_store, kinetic.dhs_evaluate(temperatures, pressures)
    )


def test_dh_not_specified_multiple_reaction():
    a = rd.Substance.from_thermo_database("ch4", "methane")
    b = rd.Substance.from_thermo_database("o2", "oxygen")
    c = rd.Substance.from_thermo_database("co2", "carbon dioxide")
    d = rd.Substance.from_thermo_database("co", "carbon monoxide")
    e = rd.Substance.from_thermo_database("h2o", "water")

    mix = rd.mix.IdealGas([a, b, c, d, e])

    def rate():
        ...

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a + 2 * b > c + 2 * e, "rate": rate},
            "r2": {"eq": a + 3 / 2 * b > d + e, "rate": rate},
            "r3": {"eq": 2 * d + b > c, "rate": rate},
        },
        kinetic_constants={},
    )

    kinetic.set_dh_function()

    cons, corr = ChemicalConstantsPackage.from_IDs(
        ["methane", "oxygen", "carbon dioxide", "carbon monoxide", "water"]
    )

    hs = cons.Hfgs
    cps = corr.HeatCapacityGases

    deltah_std1 = (hs[2] + 2 * hs[4]) - (hs[0] + 2 * hs[1])
    deltah_std2 = (hs[3] + hs[4]) - (hs[0] + 3 / 2 * hs[1])
    deltah_std3 = hs[2] - (2 * hs[3] + hs[1])

    assert kinetic.dhs_evaluate(298.15, 101325).ravel()[0] == deltah_std1
    assert kinetic.dhs_evaluate(298.15, 101325).ravel()[1] == deltah_std2
    assert kinetic.dhs_evaluate(298.15, 101325).ravel()[2] == deltah_std3

    temperatures = np.linspace(298.15, 700, 100)
    pressures = np.full(100, 101325)

    dhs_store = np.array([])

    # One by one
    for t in temperatures:
        integrals = np.array(
            [quad(cp.T_dependent_property, 298.15, t)[0] for cp in cps]
        )

        dh_corrected = np.array(
            [
                deltah_std1
                + (integrals[2] + 2 * integrals[4])
                - (integrals[0] + 2 * integrals[1]),
                deltah_std2
                + (integrals[3] + integrals[4])
                - (integrals[0] + 3 / 2 * integrals[1]),
                deltah_std3 + integrals[2] - (2 * integrals[3] + integrals[1]),
            ]
        )

        assert np.allclose(
            dh_corrected, kinetic.dhs_evaluate(t, 101325), atol=1e-10
        )

        assert (
            kinetic.dhs_evaluate(t, 101325)
            == dh_not_specified(kinetic, t, 101325)
        ).all()

        if t == 298.15:
            dhs_store = np.append(dhs_store, dh_corrected)
        else:
            dhs_store = np.vstack((dhs_store, dh_corrected))

    # vectorial call
    assert np.allclose(
        dhs_store.T, kinetic.dhs_evaluate(temperatures, pressures)
    )
