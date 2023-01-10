import numpy as np

import pytest

import reactord as rd


def test_simple_evaluation1():

    mix = rd.mix.IdealSolution(a="water", b="ethanol", c="acetone")

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

    """ethanol = rd.Substance.from_thermo_database("ethanol")
    acetic = rd.Substance.from_thermo_database("acetic acid")
    water = rd.Substance.from_thermo_database("water")
    acetate = rd.Substance.from_thermo_database("ethyl acetate")

    list_of_components = [acetic, ethanol, acetate, water]"""

    stoichiometry_single_reaction = np.array([-1, -1, 1, 1])

    def reaction_rate(concentration, temperature):
        return 10

    def reaction_rate2(concentration, temperature):
        return 15

    mix1 = rd.mix.IdealGas(
        A="acetic acid", B="ethanol", C="ethyl acetate", D="water"
    )
    mix2 = rd.mix.IdealSolution(
        A="acetic acid", B="ethanol", C="ethyl acetate", D="water"
    )

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
    for component in mix1.substances:
        raw_formation_enthalpies.append(component.formation_enthalpy_ig)
    raw_formation_enthalpies = np.array(raw_formation_enthalpies)

    # stoichiometry is rehsaped to calculate the dot product
    stoichiometry = np.reshape(kinetic1.stoichiometry, newshape=(4,))
    raw_reaction_enthalpy = np.dot(raw_formation_enthalpies, stoichiometry)

    formation_enthalpies_method = (
        kinetic1._std_reaction_enthalpies_from_formation()
    )

    assert np.allclose(raw_reaction_enthalpy, formation_enthalpies_method)


def test_standard_reaction_enthalpies():
    def reac1(concentrations, temperature):
        return 10

    def reac2(concentrations, temperature):
        return 20

    stoichiometry = np.array([[-1, -1, 1, 1], [0, -2, 0, 2]])

    a = rd.Substance.from_thermo_database("acetic acid")
    b = rd.Substance.from_thermo_database("hydrogen peroxide")
    c = rd.Substance.from_thermo_database("peracetic acid")
    d = rd.Substance.from_thermo_database("water")

    mix = rd.mix.IdealSolution(a=a, b=b, c=c, d=d)

    kinetic = rd.Kinetics(mix, [reac1, reac2], stoichiometry, "concentration")

    kinetic.std_reaction_enthalpies_init()

    enthalpy1, enthalpy2 = kinetic.std_reaction_enthalpies

    h1 = (
        c.formation_enthalpy
        + d.formation_enthalpy
        - a.formation_enthalpy
        - b.formation_enthalpy
    )

    h2 = 2 * d.formation_enthalpy - 2 * b.formation_enthalpy

    assert enthalpy1 == h1
    assert enthalpy2 == h2


def test_user_defined_reaction_enthalpies():
    def reac1(concentrations, temperature):
        return 10

    def reac2(concentrations, temperature):
        return 20

    stoichiometry = np.array([[-1, -1, 1, 1], [0, -1, 0, 1]])

    a = rd.Substance.from_thermo_database("acetic acid")
    b = rd.Substance.from_thermo_database("hydrogen peroxide")
    c = rd.Substance.from_thermo_database("peracetic acid")
    d = rd.Substance.from_thermo_database("water")

    mix = rd.mix.IdealSolution(a=a, b=b, c=c, d=d)

    kinetic = rd.Kinetics(
        mix,
        [reac1, reac2],
        stoichiometry,
        "concentration",
        reaction_enthalpies=[1000, 2000],
    )

    kinetic.std_reaction_enthalpies_init()

    assert kinetic.std_reaction_enthalpies is None

    assert kinetic.reaction_enthalpies(500, 101325) == [1000, 2000]


def test_enthalpies_correction():
    def reac1(concentrations, temperature):
        return 10

    def reac2(concentrations, temperature):
        return 20

    stoichiometry = np.array([[-1, -1, 1, 1], [0, -2, 0, 2]])

    a = rd.Substance.from_thermo_database("acetic acid")
    b = rd.Substance.from_thermo_database("hydrogen peroxide")
    c = rd.Substance.from_thermo_database("peracetic acid")
    d = rd.Substance.from_thermo_database("water")

    mix = rd.mix.IdealSolution(a=a, b=b, c=c, d=d)

    kinetic = rd.Kinetics(mix, [reac1, reac2], stoichiometry, "concentration")

    kinetic.std_reaction_enthalpies_init()

    enthalpy1, enthalpy2 = kinetic.reaction_enthalpies(330, 101325)

    h1 = (
        c.formation_enthalpy
        + c.heat_capacity_liquid_dt_integral(298.15, 330)
        + d.formation_enthalpy
        + d.heat_capacity_liquid_dt_integral(298.15, 330)
        - (
            a.formation_enthalpy
            + a.heat_capacity_liquid_dt_integral(298.15, 330)
        )
        - (
            b.formation_enthalpy
            + b.heat_capacity_liquid_dt_integral(298.15, 330)
        )
    )

    h2 = 2 * (
        d.formation_enthalpy + d.heat_capacity_liquid_dt_integral(298.15, 330)
    ) - 2 * (
        b.formation_enthalpy + b.heat_capacity_liquid_dt_integral(298.15, 330)
    )

    assert np.allclose(enthalpy1, h1)
    assert np.allclose(enthalpy2, h2)


def test_making_it_explode():
    def reac1(concentrations, temperature):
        return 10

    def reac2(concentrations, temperature):
        return 20

    stoichiometry = np.array([[-1, -1, 1, 1], [0, -2, 0, 2]])

    a = rd.Substance.from_thermo_database("acetic acid")
    b = rd.Substance.from_thermo_database("hydrogen peroxide")
    c = rd.Substance.from_thermo_database("peracetic acid")
    d = rd.Substance.from_thermo_database("water")

    mix = rd.mix.IdealSolution(a=a, b=b, c=c, d=d)

    with pytest.raises(TypeError):
        kinetic = rd.Kinetics(
            "Muchachos", [reac1, reac2], stoichiometry, "concentration"
        )

    with pytest.raises(TypeError):
        kinetic = rd.Kinetics(
            3.14, [reac1, reac2], stoichiometry, "concentration"
        )

    mix2 = rd.mix.IdealSolution(a="methane")

    with pytest.raises(IndexError):
        kinetic = rd.Kinetics(
            mix2, [reac1, reac2], stoichiometry, "concentration"
        )

    with pytest.raises(NotImplementedError):
        kinetic = rd.Kinetics(
            mix, [reac1, reac2], stoichiometry, "concentration"
        )

        kinetic.std_reaction_enthalpies = [10, 20, 30, 0]
