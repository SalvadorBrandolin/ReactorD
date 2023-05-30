import pytest

import reactord as rd
import reactord.flowreactors.stationary_1d.pfr as pfr


def test_raises():
    def r_rate():
        ...

    a = rd.Substance("A")
    b = rd.Substance("B")
    c = rd.Substance("C")

    mix = rd.mix.IdealGas([a, b, c])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c, "rate": r_rate}},
        kinetic_constants={},
    )

    mb = pfr.mass_balances.MolarFlow({"A": 10, "B": 10})
    eb = pfr.energy_balances.Isothermic(303.15)
    pb = pfr.pressure_balances.Isobaric(101325)

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    with pytest.raises(ValueError):
        reactor.initial_profile_builder()

    mb2 = pfr.mass_balances.MolarFlow({"A": 10, "B": 10, "C": 20}, {"A": 2})

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb2,
        energy_balance=eb,
        pressure_balance=pb,
    )

    with pytest.raises(ValueError):
        reactor.initial_profile_builder()
