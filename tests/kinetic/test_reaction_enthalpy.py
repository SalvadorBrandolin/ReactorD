import numpy as np

import reactord as rd
from reactord.kinetic.reaction_enthalpy import dh_not_specified, dh_specified

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

    kinetic._init_dh_function()

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

    kinetic.init_dh_function()

    cons, corr = ChemicalConstantsPackage.from_IDs(["methane", "oxygen", "carbon dioxide", "water"])
    
    hs = cons.Hfgs
    deltah = hs
