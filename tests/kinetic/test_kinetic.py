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
        return a * np.exp(-e / R / t) * c["etol"] * c["ach"]

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c + d, "rate": r1_rate}},
        kinetic_constants={"a": 10, "e": 10000},
        rates_argument="concentration",
    )

    temperatures = np.linspace(300, 600, 50)
    pressures = np.linspace(100_000, 300_000, 50)
    moles = np.array([np.linspace(0.1, 1, 50), np.linspace(0, 4, 50), np.linspace(0, 0.8, 50), np.linspace(0, 0.3, 50)])
    
    # One by one
    for m, t, p in zip(moles.T, temperatures, pressures):
        import ipdb
        ipdb.set_trace()
        z_rd = kinetic.mix.mole_fractions(m)
        r_reactord = kinetic.evaluate(z_rd, t, p)
        
        # by hand
        z_bh = m / np.sum(m)
        rho = p / R / t
        c = z_bh * rho
        r_byhand = 10 * np.exp(-10000 / R / t) * c[0] * c[1]
        
        assert r_reactord == r_byhand
