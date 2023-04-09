from IPython.display import Math

from typing import Union

import numpy as np

import reactord as rd

from sympy import Eq, Equality, Symbol, latex, symbols

import pytest


def test_errors():
    a = rd.Substance("a")
    b = rd.Substance("b")

    # name error
    with pytest.raises(TypeError):
        rd.Substance(name=2001)
    with pytest.raises(TypeError):
        rd.Substance(name=["methane"])

    # __mul__
    with pytest.raises(TypeError):
        a * b
    with pytest.raises(TypeError):
        a * "-team"

    # __sum__
    with pytest.raises(TypeError):
        a + 3
    with pytest.raises(TypeError):
        b + "atistuta"

    # __gt__
    with pytest.raises(TypeError):
        a + b > 3
    with pytest.raises(TypeError):
        b + a > "nana"


def test_initial_symbol():
    a = rd.Substance("a")
    b = rd.Substance("b")
    c = rd.Substance("methane")
    
    # Initial Symbolic
    assert a._expression == symbols("a")
    assert a._chem_equality is None
    assert a._names == np.array(["a"])

    assert b._expression == symbols("b")
    assert b._chem_equality is None
    assert b._names == np.array(["b"])
    
    assert c._expression == symbols("methane")
    assert c._chem_equality is None
    assert c._names == np.array(["methane"])

@pytest.mark.parametrize("coeff", [1, 2, 1/2, 5/3, 0.658])
def test_mul_operand(coeff: Union[int, float]):
    a = rd.Substance("a")
        
    expr1 = coeff * a
    assert expr1._expression == coeff * symbols(f"a")
    assert expr1._chem_equality is None
    assert expr1._names == np.array(["a"])
    
    expr2 = a * coeff
    assert expr2._expression == coeff * symbols(f"a")
    assert expr2._chem_equality is None
    assert expr2._names == np.array(["a"])
    
    expr3 = coeff * a * coeff
    assert expr3._expression == coeff**2 * symbols(f"a")
    assert expr3._chem_equality is None
    assert expr3._names == np.array(["a"])

def test_add_opereand1():
    a = rd.Substance("a")
    b = rd.Substance("b")
    
    a_sym = symbols("a")
    b_sym = symbols("b")
    
    expr1 = a + b
    assert expr1._expression == a_sym + b_sym
    assert expr1._chem_equality is None
    assert (expr1._names == np.array(["a", "b"])).all()
    
    expr2 = a + 2 * b
    assert expr2._expression == a_sym + 2 * b_sym
    assert expr2._chem_equality is None
    assert (expr2._names == np.array(["a", "b"])).all()
    
    expr3 = a + 2 * b * 5 / 2
    assert expr3._expression == a_sym + 5 * b_sym
    assert expr3._chem_equality is None
    assert (expr3._names == np.array(["a", "b"])).all()
    
    expr4 = 2 * a + b
    assert expr4._expression == 2 * a_sym + b_sym
    assert expr4._chem_equality is None
    assert (expr4._names == np.array(["a", "b"])).all()
    
    expr5 = 2 * a * 3 / 5 + b
    assert expr5._expression == 6 / 5 * a_sym + b_sym
    assert expr5._chem_equality is None
    assert (expr5._names == np.array(["a", "b"])).all()
    
    expr6 = 2 * a * 3 + 5 * b
    assert expr6._expression == 6 * a_sym + 5 * b_sym
    assert expr6._chem_equality is None
    assert (expr6._names == np.array(["a", "b"])).all()
    
    expr7 = 2 * a * 3 + 5 * b * 2
    assert expr7._expression == 6 * a_sym + 10 * b_sym
    assert expr7._chem_equality is None
    assert (expr7._names == np.array(["a", "b"])).all()
    