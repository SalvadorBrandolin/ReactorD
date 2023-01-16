import numpy as np

import reactord as rd


# Example: methane combustion reaction
# CH4 + 2 O2 ----> 2 H2O + CO2
# A  + 2B  ----> 2C + D


# We create the substance objects:
A = rd.Substance.from_thermo_database("methane")
B = rd.Substance.from_thermo_database("oxygen")
C = rd.Substance.from_thermo_database("water")
D = rd.Substance.from_thermo_database("co2")

# Then we create the mixture:
mixture = rd.mix.IdealGas(A=A, B=B, C=C, D=D)

# and the stoichiometry matrix:
stoichiometry = np.array([[-1, -2, 2, 1]])

# A function for the reaction rate is defined:


def rate1(concentrations, temperature):
    return 10000 * concentrations[0]


# Kinetics is instantiated and "concentration" are being used to make
# the calculations
kinetics = rd.Kinetics(mixture, [rate1], stoichiometry, "concentration")

# An instance of the stationary plug flow reactor class is instantiated
# at isothermic and isobaric conditions:
PFR = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
    mix=mixture,
    list_of_reactions=[rate1],
    stoichiometry=stoichiometry,
    kinetic_argument="partial_pressure",
    reactor_dim_minmax=[0, 1],
    transversal_area=1,
    isothermic_temperature=350,
    isobaric_pressure=5 * 101325,
    molar_flow_in={
        "methane": 100,
        "oxygen": 100,
        "water": 0,
        "carbon dioxide": 0,
    },
)

# Simulation:
solution = PFR.simulate()
