import pytest

import reactord as rd


def test_repr():
    pb = rd.flowreactors.stationary_1d.pfr.pressure_balances.Ergun(
        {"in": 101325}, porosity=0.1, particle_diameter=1
    )

    text = (
        r"\frac{dP}{dz}=-\frac{G}{{\rho}D_p}\left(\frac{1-\phi}{\phi^3}"
        r"\right)\left[\frac{150(1-\phi)\mu}{D_p}+1.75G\right]"
    )

    assert pb.__repr__() == text


def test_raises():
    with pytest.raises(ValueError):
        rd.flowreactors.stationary_1d.pfr.pressure_balances.Ergun(
            pressure={"in": 300_000, "out": 100_000},
            porosity=0.4,
            particle_diameter=1,
        )
