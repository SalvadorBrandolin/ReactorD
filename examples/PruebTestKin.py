import numpy as np

import pytest

import reactord as rd


def test_user_reaction_enthalpies():
    ethanol = rd.Substance.from_thermo_database("ethanol")
    acetic = rd.Substance.from_thermo_database("acetic acid")
    water = rd.Substance.from_thermo_database("water")
    acetate = rd.Substance.from_thermo_database("ethyl acetate")

    list_of_components = [acetic, ethanol, acetate, water]
    stoichiometry = np.array([-1, -1, 1, 1])

    def reaction_rate(concentration, temperature):
        return 10

    mix2 = rd.mix.IdealSolution(list_of_components)
    wrong_reaction_enthalpies = np.array([10, 10])

    with pytest.raises(IndexError):
        rd.Kinetics(
            mix=mix2,
            list_of_reactions=[reaction_rate],
            stoichiometry=stoichiometry,
            kinetic_argument="concentration",
            reaction_enthalpies=wrong_reaction_enthalpies,
        )

    with pytest.raises(IndexError):
        rd.Kinetics(
            mix=mix2,
            list_of_reactions=[reaction_rate],
            stoichiometry=stoichiometry,
            kinetic_argument="concentration",
            reaction_enthalpies=wrong_reaction_enthalpies,
        )
