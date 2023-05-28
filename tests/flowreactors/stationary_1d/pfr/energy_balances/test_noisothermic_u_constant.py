import numpy as np

import pytest

from reactord.flowreactors.stationary_1d.pfr.energy_balances import (
    NoIsothermicUConstant,
)


def test_raises():
    with pytest.raises(NotImplementedError):
        NoIsothermicUConstant(
            temperature_in_or_out={"in": 200, "out": 400},
            refrigerant_in_temperature=800,
            heat_exchange_coefficient=110,
            refrigerant_mix=None,
            refrigerant_composition=np.array([1]),
            refrigerant_pressure=101325,
            refrigerant_molar_flow=0.10,
            heat_exchange_arrangement="co-current",
        )

    with pytest.raises(ValueError):
        NoIsothermicUConstant(
            temperature_in_or_out={"in": 200},
            refrigerant_in_temperature=800,
            heat_exchange_coefficient=110,
            refrigerant_mix=None,
            refrigerant_composition=np.array([1]),
            refrigerant_pressure=101325,
            refrigerant_molar_flow=0.10,
            heat_exchange_arrangement="contracorriente",
        )


def test_repr():
    eb1 = NoIsothermicUConstant(
        temperature_in_or_out={"in": 200},
        refrigerant_in_temperature=800,
        heat_exchange_coefficient=110,
        refrigerant_mix=None,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.10,
        heat_exchange_arrangement="co-current",
    )

    latex1 = (
        r"\frac{1}{a_t}\frac{dT}{dz}=\frac{Ua(T_a-T)+\sum\Delta{H_j}r_{j}}"
        r"{{c_p}_{mix}{\sum}F_i}"
    )

    latex2 = (
        r"\frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T-T_r)}{{c_p}_{r}{\sum"
        r"}F_i}}"
    )

    latex3 = (
        r"\frac{1}{a_t}\frac{dT_r}{dz}=\frac{Ua(T_r-T)}{{c_p}_{r}{\sum"
        r"}F_i}}"
    )

    assert (latex1, latex2) == eb1.__repr__()

    eb2 = NoIsothermicUConstant(
        temperature_in_or_out={"in": 200},
        refrigerant_in_temperature=800,
        heat_exchange_coefficient=110,
        refrigerant_mix=None,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.10,
        heat_exchange_arrangement="counter-current",
    )

    assert (latex1, latex3) == eb2.__repr__()
