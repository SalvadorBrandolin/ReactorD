import pytest

import reactord.flowreactors.stationary_1d.pfr as pfr


def test_raises():
    with pytest.raises(NotImplementedError):
        pfr.energy_balances.NoIsothermicAllConstant(
            temperature_in_or_out={"in": 300, "out": 400},
            refrigerant_in_temperature=500,
            heat_exchange_coefficient=110,
        )


def test_repr():
    eb = pfr.energy_balances.NoIsothermicAllConstant(
        temperature_in_or_out={"in": 300},
        refrigerant_in_temperature=500,
        heat_exchange_coefficient=110,
    )

    latex1 = (
        r"\frac{1}{a_t}\frac{dT}{dz}=\frac{Ua(T_a-T)+\sum\Delta{H_j}r_{j}}{{c_"
        r"p}_{mix}{\sum}F_i}"
    )
    latex2 = r"\frac{dT_r}{dz}=0"

    text = (latex1, latex2)

    assert eb.__repr__() == text
