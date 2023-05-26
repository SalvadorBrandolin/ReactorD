import numpy as np

import pytest

import reactord as rd

from scipy.constants import R


def test_repr():
    def rate():
        ...

    a = rd.Substance("A")
    b = rd.Substance("B")
    c = rd.Substance("C")
    d = rd.Substance("D")
    e = rd.Substance("E")

    mix = rd.mix.IdealGas([a, b, c, d, e])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a > b, "rate": rate},
            "r2": {"eq": b + 2 * c > d, "rate": rate},
            "r3": {"eq": 2 * d + a > 2 * e, "rate": rate},
        },
        kinetic_constants={},
    )

    text = (
        "Mixture's substances: \n"
        "  * A \n"
        "  * B \n"
        "  * C \n"
        "  * D \n"
        "  * E \n"
        "\n"
        "System's reactions: \n"
        "r1: A \\rightarrow B \n"
        "r2: B + 2 C \\rightarrow D \n"
        "r3: A + 2 D \\rightarrow 2 E \n"
    )

    assert text == kinetic.__repr__()


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
    moles = np.array(
        [
            np.linspace(0.1, 1, 50),
            np.linspace(0, 4, 50),
            np.linspace(0, 0.8, 50),
            np.linspace(0, 0.3, 50),
        ]
    )

    r_store = np.array([])

    # One by one
    for m, t, p in zip(moles.T, temperatures, pressures):
        z_rd = kinetic.mix.mole_fractions(m)
        r_reactord = kinetic.evaluate(z_rd, t, p)

        # by hand
        z_bh = m / np.sum(m)
        rho = p / R / t
        c = z_bh * rho
        r_byhand = 10 * np.exp(-10000 / R / t) * c[0] * c[1]

        assert np.allclose(r_reactord, r_byhand, atol=1e-10)

        r_store = np.append(r_store, r_byhand)

    # Vectorial call
    zs = kinetic.mix.mole_fractions(moles)
    assert np.allclose(
        r_store, kinetic.evaluate(zs, temperatures, pressures), atol=1e-10
    )


def test_one_reaction2():
    """NOT REAL KINETIC PARAMETERS."""
    a = rd.Substance.from_thermo_database("ch4", "methane")
    b = rd.Substance.from_thermo_database("h2o", "water")
    c = rd.Substance.from_thermo_database("co", "carbon monoxide")
    d = rd.Substance.from_thermo_database("h2", "hydrogen")

    mix = rd.mix.IdealGas([a, b, c, d])

    def r1_rate(c, t, cons):
        a_d, e_d, a_i, e_i = cons["a_d"], cons["e_d"], cons["a_i"], cons["e_i"]
        rate = (
            a_d * np.exp(-e_d / R / t) * c["ch4"] * c["h2o"]
            - a_i * np.exp(-e_i / R / t) * c["co"] * c["h2"] ** 3
        )
        return rate

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c + 3 * d, "rate": r1_rate}},
        kinetic_constants={"a_d": 2, "e_d": 8000, "a_i": 1, "e_i": 10000},
        rates_argument="partial pressure",
    )

    temperatures = np.linspace(300, 600, 50)
    pressures = np.linspace(100_000, 300_000, 50)
    moles = np.array(
        [
            np.linspace(0.1, 1, 50),
            np.linspace(0, 4, 50),
            np.linspace(0, 0.8, 50),
            np.linspace(0, 0.3, 50),
        ]
    )

    r_store = np.array([])

    # One by one
    for m, t, p in zip(moles.T, temperatures, pressures):
        z_rd = kinetic.mix.mole_fractions(m)
        r_reactord = kinetic.evaluate(z_rd, t, p)

        # by hand
        z_bh = m / np.sum(m)
        pp = z_bh * p
        r_byhand = (
            2 * np.exp(-8000 / R / t) * pp[0] * pp[1]
            - 1 * np.exp(-10000 / R / t) * pp[2] * pp[3] ** 3
        )

        assert np.allclose(r_reactord, r_byhand, atol=1e-10)

        r_store = np.append(r_store, r_byhand)

    # Vectorial call
    zs = kinetic.mix.mole_fractions(moles)
    assert np.allclose(
        r_store, kinetic.evaluate(zs, temperatures, pressures), atol=1e-10
    )


def test_multiple_reactions():
    """NOT REAL KINETIC PARAMETERS."""
    ch4 = rd.Substance.from_thermo_database("ch4", "methane")
    c2h6 = rd.Substance.from_thermo_database("c2h6", "ethane")
    h2o = rd.Substance.from_thermo_database("h2o", "water")
    co = rd.Substance.from_thermo_database("co", "carbon monoxide")
    h2 = rd.Substance.from_thermo_database("h2", "hydrogen")
    co2 = rd.Substance.from_thermo_database("co2", "carbon dioxide")
    n2 = rd.Substance.from_thermo_database("n2", "nitrogen")

    mix = rd.mix.IdealGas([ch4, c2h6, h2o, co, h2, co2, n2])

    def rate1(c, t, cons):
        ad1, ed1, ai1, ei1 = cons["ad1"], cons["ed1"], cons["ai1"], cons["ei1"]
        rate = (
            ad1 * np.exp(-ed1 / R / t) * c["ch4"] * c["h2o"]
            - ai1 * np.exp(-ei1 / R / t) * c["co"] * c["h2"] ** 3
        )
        return rate

    def rate2(c, t, cons):
        ad2, ed2, ai2, ei2 = cons["ad2"], cons["ed2"], cons["ai2"], cons["ei2"]
        rate = (
            ad2 * np.exp(-ed2 / R / t) * c["co"] * c["h2o"]
            - ai2 * np.exp(-ei2 / R / t) * c["co2"] * c["h2"]
        )
        return rate

    def rate3(c, t, cons):
        ad3, ed3, ai3, ei3 = cons["ad3"], cons["ed3"], cons["ai3"], cons["ei3"]
        rate = (
            ad3 * np.exp(-ed3 / R / t) * c["c2h6"] * c["co2"] ** 2
            - ai3 * np.exp(-ei3 / R / t) * c["co"] ** 4 * c["h2"] ** 3
        )
        return rate

    kinetic_concentration = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": ch4 + h2o > co + 3 * h2, "rate": rate1},
            "r2": {"eq": co + h2o > co2 + h2, "rate": rate2},
            "r3": {"eq": c2h6 + 2 * co2 > 4 * co + 3 * h2, "rate": rate3},
        },
        kinetic_constants={
            "ad1": 1,
            "ed1": 5000,
            "ai1": 0.5,
            "ei1": 10000,
            "ad2": 2,
            "ed2": 6000,
            "ai2": 1.2,
            "ei2": 8000,
            "ad3": 3,
            "ed3": 3000,
            "ai3": 2,
            "ei3": 4500,
        },
        rates_argument="concentration",
    )

    kinetic_partial_pressure = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": ch4 + h2o > co + 3 * h2, "rate": rate1},
            "r2": {"eq": co + h2o > co2 + h2, "rate": rate2},
            "r3": {"eq": c2h6 + 2 * co2 > 4 * co + 3 * h2, "rate": rate3},
        },
        kinetic_constants={
            "ad1": 1,
            "ed1": 5000,
            "ai1": 0.5,
            "ei1": 10000,
            "ad2": 2,
            "ed2": 6000,
            "ai2": 1.2,
            "ei2": 8000,
            "ad3": 3,
            "ed3": 3000,
            "ai3": 2,
            "ei3": 4500,
        },
        rates_argument="partial pressure",
    )

    temperatures = np.linspace(300, 600, 1000)
    pressures = np.linspace(100_000, 300_000, 1000)
    moles = np.array(
        [
            np.linspace(0.1, 1, 1000),
            np.linspace(0, 4, 1000),
            np.linspace(0, 0.8, 1000),
            np.linspace(0, 0.3, 1000),
            np.linspace(0.2, 0.7, 1000),
            np.linspace(0.5, 1, 1000),
            np.full(1000, 5),
        ]
    )

    mix = rd.mix.IdealGas([ch4, c2h6, h2o, co, h2, co2, n2])

    # =========================================================================
    # Concentration
    # =========================================================================
    r_store = np.array([])

    # One by one
    for m, t, p in zip(moles.T, temperatures, pressures):
        z_rd = kinetic_concentration.mix.mole_fractions(m)
        r_reactord = kinetic_concentration.evaluate(z_rd, t, p)

        # by hand
        z_bh = m / np.sum(m)
        rho = p / R / t
        c = z_bh * rho
        r_byhand = np.array(
            [
                1 * np.exp(-5000 / R / t) * c[0] * c[2]
                - 0.5 * np.exp(-10000 / R / t) * c[3] * c[4] ** 3,
                2 * np.exp(-6000 / R / t) * c[3] * c[2]
                - 1.2 * np.exp(-8000 / R / t) * c[5] * c[4],
                3 * np.exp(-3000 / R / t) * c[1] * c[5] ** 2
                - 2 * np.exp(-4500 / R / t) * c[3] ** 4 * c[4] ** 3,
            ]
        )

        assert np.allclose(r_reactord, r_byhand, atol=1e-10)

        if t == 300:
            r_store = r_byhand
        else:
            r_store = np.vstack((r_store, r_byhand))

    # Vectorial call
    zs = kinetic_concentration.mix.mole_fractions(moles)
    assert np.allclose(
        r_store.T, kinetic_concentration.evaluate(zs, temperatures, pressures)
    )

    # =========================================================================
    # Partial pressure
    # =========================================================================
    r_store = np.array([])

    # One by one
    for m, t, p in zip(moles.T, temperatures, pressures):
        z_rd = kinetic_partial_pressure.mix.mole_fractions(m)
        r_reactord = kinetic_partial_pressure.evaluate(z_rd, t, p)

        # by hand
        z_bh = m / np.sum(m)
        pp = z_bh * p
        r_byhand = np.array(
            [
                1 * np.exp(-5000 / R / t) * pp[0] * pp[2]
                - 0.5 * np.exp(-10000 / R / t) * pp[3] * pp[4] ** 3,
                2 * np.exp(-6000 / R / t) * pp[3] * pp[2]
                - 1.2 * np.exp(-8000 / R / t) * pp[5] * pp[4],
                3 * np.exp(-3000 / R / t) * pp[1] * pp[5] ** 2
                - 2 * np.exp(-4500 / R / t) * pp[3] ** 4 * pp[4] ** 3,
            ]
        )

        assert np.allclose(r_reactord, r_byhand, atol=1e-10)

        if t == 300:
            r_store = r_byhand
        else:
            r_store = np.vstack((r_store, r_byhand))

    # Vectorial call
    zs = kinetic_concentration.mix.mole_fractions(moles)
    assert np.allclose(
        r_store.T,
        kinetic_partial_pressure.evaluate(zs, temperatures, pressures),
    )


def test_raises():
    a = rd.Substance("a")
    b = rd.Substance("b")
    c = rd.Substance("c")

    mix = rd.mix.IdealGas([a, b, c])

    def rate():
        ...

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a + b > c, "rate": rate, "DH": 1},
            "r2": {"eq": 2 * a + b > 3 * c, "rate": rate, "DH": 2},
        },
        kinetic_constants={},
    )

    assert (kinetic.user_r_dh == np.array([1, 2])).all()

    with pytest.raises(ValueError):
        kinetic.user_r_dh = np.array([3, 4])

    with pytest.raises(NotImplementedError):
        rd.Kinetic(
            mix=mix,
            reactions={
                "r1": {"eq": a + b > c, "rate": rate},
                "r2": {"eq": 2 * a + b > 3 * c, "rate": rate, "DH": 2},
            },
            kinetic_constants={},
        )
