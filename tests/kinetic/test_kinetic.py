import numpy as np

import reactord as rd

from scipy.constants import R


def test_one_reaction1():
    """NOT REAL KINETIC PARAMETERS."""

    a = rd.Substance.from_thermo_database("etol", "ethanol")
    b = rd.Substance.from_thermo_database("ach", "acetic acid")
    c = rd.Substance.from_thermo_database("etac", "ethyl acetate")
    d = rd.Substance.from_thermo_database("wt", "water")

    mix = rd.mix.IdealGas([a, b, c, d])

    def r1_rate(c, t, cons):
        a, e = cons["a"], cons["e"]
        return a * np.exp(-e / 8.314 / t) * c["etol"] * c["ach"]

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c + d, "rate": r1_rate}},
        kinetic_constants={"a": 10, "e": 10000},
        rates_argument="concentration",
    )

    assert (kinetic.stoichiometry == np.array([-1, -1, 1, 1])).all()
    assert np.ndim(kinetic.stoichiometry) == 2
    assert np.shape(kinetic.stoichiometry) == (1, 4)

    p = 101325
    t = 303.15
    x = [0.25, 0.25, 0.25, 0.25]

    density = p / R / t
