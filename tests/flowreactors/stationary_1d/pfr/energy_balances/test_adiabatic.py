import pytest

import reactord.flowreactors.stationary_1d.pfr as pfr


def test_raises():
    with pytest.raises(NotImplementedError):
        pfr.energy_balances.Adiabatic({"in": 303.15, "out": 600})
