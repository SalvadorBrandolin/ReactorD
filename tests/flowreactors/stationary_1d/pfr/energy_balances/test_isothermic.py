import reactord.flowreactors.stationary_1d.pfr as pfr


def test_repr():
    eb = pfr.energy_balances.Isothermic(298.15)

    text = (r"\frac{dT}{dz}=0", r"\frac{dT_r}{dz}=0")

    assert text == eb.__repr__()
