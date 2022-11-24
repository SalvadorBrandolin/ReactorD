import numpy as np

import pytest

import reactord as rd


def test_simple_evaluation1():

    substance_a = rd.Substance.from_thermo_database("water")
    substance_b = rd.Substance.from_thermo_database("ethanol")
    substance_c = rd.Substance.from_thermo_database("acetone")

    mix = rd.mix.IdealSolution([substance_a, substance_b, substance_c])

    def law1(concentrations, temperature):
        return 10

    def law2(concentrations, temperature):
        return 8

    stoichiometry = np.array([[-1, 1, 0], [0, -1, 1]])

    kinetics = rd.Kinetics(
        mix=mix,
        list_of_reactions=[law1, law2],
        stoichiometry=stoichiometry,
        kinetic_argument="concentration",
    )

    compositions = np.array(
        [[1, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.33, 0.46, 0.35]]
    )

    temperature = np.array([300, 400, 500, 600])

    pressure = 101325

    for z in compositions:
        for t in temperature:
            c_subs_r_rates, c_r_rates = kinetics.kinetic_eval(z, t, pressure)

            assert c_subs_r_rates[0] == -10.0
            assert c_subs_r_rates[1] == 2.0
            assert c_subs_r_rates[2] == 8.0

            assert c_r_rates[0] == 10.0
            assert c_r_rates[1] == 8.0


def test_user_reaction_enthalpies():

    ethanol = rd.Substance.from_thermo_database("ethanol")
    acetic = rd.Substance.from_thermo_database("acetic acid")
    water = rd.Substance.from_thermo_database("water")
    acetate = rd.Substance.from_thermo_database("ethyl acetate")

    list_of_components = [acetic, ethanol, acetate, water]

    stoichiometry_single_reaction = np.array([-1, -1, 1, 1])

    def reaction_rate(concentration, temperature):
        return 10

    def reaction_rate2(concentration, temperature):
        return 15

    mix1 = rd.mix.IdealGas(list_of_components)
    mix2 = rd.mix.IdealSolution(list_of_components)

    kinetic1 = rd.Kinetics(
        mix1,
        list_of_reactions=[reaction_rate],
        stoichiometry=stoichiometry_single_reaction,
        kinetic_argument="concentration",
        reaction_enthalpies=np.array([5]),
    )

    kinetic2 = rd.Kinetics(
        mix=mix2,
        list_of_reactions=[reaction_rate],
        stoichiometry=stoichiometry_single_reaction,
        kinetic_argument="partial_pressure",
        reaction_enthalpies=np.array([5]),
    )

    assert kinetic1.reaction_enthalpies(1250, 2) == np.array([5])
    assert kinetic2.reaction_enthalpies(1250, 2) == np.array([5])

    wrong_reaction_enthalpies = np.array([10, 10])

    # Test there is one reaction enthalpy per reaction
    with pytest.raises(IndexError):
        rd.Kinetics(
            mix=mix2,
            list_of_reactions=[reaction_rate],
            stoichiometry=stoichiometry_single_reaction,
            kinetic_argument="concentration",
            reaction_enthalpies=wrong_reaction_enthalpies,
        )

    # Test whether there is one reaction in list_of_reactions
    # per row in stoichiometry
    correct_reaction_enthalpies = [20, 20]
    with pytest.raises(IndexError):
        rd.Kinetics(
            mix=mix2,
            list_of_reactions=[reaction_rate, reaction_rate2],
            stoichiometry=stoichiometry_single_reaction,
            kinetic_argument="partial_pressure",
            reaction_enthalpies=correct_reaction_enthalpies,
        )

    with pytest.raises(ValueError):
        rd.Kinetics(
            mix=mix2,
            list_of_reactions=[reaction_rate],
            stoichiometry=stoichiometry_single_reaction,
            kinetic_argument="Invalid_Argument",
            reaction_enthalpies=([25]),
        )

    # Test for _std_reaction_enthalpies_from_formation method
    raw_formation_enthalpies = []

    # Individual formation enthalpies are retrieved and the list
    # is converted to a numpy array
    for component in list_of_components:
        raw_formation_enthalpies.append(component.formation_enthalpy_ig)
    raw_formation_enthalpies = np.array(raw_formation_enthalpies)

    # stoichiometry is rehsaped to calculate the dot product
    stoichiometry = np.reshape(kinetic1.stoichiometry, newshape=(4,))
    raw_reaction_enthalpy = np.dot(raw_formation_enthalpies, stoichiometry)

    formation_enthalpies_method = (
        kinetic1._std_reaction_enthalpies_from_formation()
    )

    assert np.allclose(raw_reaction_enthalpy, formation_enthalpies_method)
