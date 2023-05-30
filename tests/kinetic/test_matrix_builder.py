import numpy as np

import reactord as rd


def test_one_reaction1():
    a = rd.Substance("a")
    b = rd.Substance("b")
    c = rd.Substance("c")
    d = rd.Substance("d")

    mix = rd.mix.IdealGas([a, b, c, d])

    def r1_rate(c, t, cons):
        ...

    # =========================================================================

    kinetic1 = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c + d, "rate": r1_rate}},
        kinetic_constants={},
        rates_argument="concentration",
    )
    assert (kinetic1.stoichiometry == np.array([-1, -1, 1, 1])).all()
    assert np.ndim(kinetic1.stoichiometry) == 2
    assert np.shape(kinetic1.stoichiometry) == (1, 4)

    # =========================================================================

    kinetic2 = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": 2 * a + 3 * b > 4 * c + 5 * d, "rate": r1_rate}
        },
        kinetic_constants={},
        rates_argument="concentration",
    )
    assert (kinetic2.stoichiometry == np.array([-2, -3, 4, 5])).all()
    assert np.ndim(kinetic2.stoichiometry) == 2
    assert np.shape(kinetic2.stoichiometry) == (1, 4)

    # =========================================================================

    kinetic3 = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a / 2 + 4 * b > 2 * c + 1 * d, "rate": r1_rate}
        },
        kinetic_constants={},
        rates_argument="concentration",
    )
    assert (kinetic3.stoichiometry == np.array([-0.5, -4, 2, 1])).all()
    assert np.ndim(kinetic3.stoichiometry) == 2
    assert np.shape(kinetic3.stoichiometry) == (1, 4)


def test_multiple_reaction():
    a = rd.Substance("a")
    b = rd.Substance("b")
    c = rd.Substance("c")
    d = rd.Substance("d")
    e = rd.Substance("e")
    f = rd.Substance("f")
    g = rd.Substance("g")

    mix = rd.mix.IdealGas([a, b, c, d, e, f, g])

    def rate():
        ...

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={
            "r1": {"eq": a + b > c, "rate": rate},
            "r2": {"eq": 2 * c + d / 4 > 1 / 5 * e, "rate": rate},
            "r3": {"eq": 5 * e + f > 2 * a + b, "rate": rate},
        },
        kinetic_constants={},
        rates_argument="concentration",
    )

    manual_stoichiometry = np.array(
        [
            [-1, -1, 1, 0, 0, 0, 0],
            [0, 0, -2, -1 / 4, 1 / 5, 0, 0],
            [2, 1, 0, 0, -5, -1, 0],
        ]
    )

    assert np.shape(kinetic.stoichiometry) == (3, 7)
    assert (kinetic.stoichiometry == manual_stoichiometry).all()
