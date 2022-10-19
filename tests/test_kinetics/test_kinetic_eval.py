import numpy as np

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
