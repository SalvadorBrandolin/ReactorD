import numpy as np

import pytest

import reactord as rd


def test_argument1():
    names = ["nm1", "nm2", "nm3", "nm4"]
    argument = rd.kinetic.CompositionalArgument(names)
    argument.values = np.array([0, 1, 2, 3])

    for idx, name in enumerate(names):
        assert argument[name] == idx


def test_argument2():
    names = ["name1", "name2", "name3"]

    values = np.zeros((3, 50))
    values[1, :] = 1
    values[2, :] = 2

    argument = rd.kinetic.CompositionalArgument(names)
    argument.values = values

    for idx, name in enumerate(names):
        assert (argument[name] == np.full(50, idx)).all()


def test_error():
    names = ["methane", "water", "carbon_dioxide"]

    argument = rd.kinetic.CompositionalArgument(names)
    argument.values = np.zeros((3, 10))

    with pytest.raises(KeyError):
        argument["hydrogen"]
