import reactord as rd


def test_repr():
    pb = rd.flowreactors.stationary_1d.pfr.pressure_balances.Isobaric(101325)

    assert pb.__repr__() == r"\frac{dP}{dz}=0"
