import numpy as np

import reactord as rd

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


# mix1: Ideal Gas
raw_formation_enthalpies = []

for component in list_of_components:
    raw_formation_enthalpies.append(component.formation_enthalpy_ig)

raw_formation_enthalpies = np.array(raw_formation_enthalpies)


print(
    raw_formation_enthalpies,
    type(raw_formation_enthalpies),
    np.shape(raw_formation_enthalpies),
    np.size(raw_formation_enthalpies),
)
print(
    kinetic1.stoichiometry,
    type(kinetic1.stoichiometry),
    np.shape(kinetic1.stoichiometry),
    np.size(kinetic1.stoichiometry),
    "\n",
)

esteq = np.reshape(kinetic1.stoichiometry, newshape=(4,))
# Z = np.reshape(raw_formation_enthalpies, newshape=(4,))
# print ("Antes de reshape: ", raw_formation_enthalpies)
# print ("Despues de reshape: ", Z, type(Z))
# print("Despues de reshape: ", esteq, type(esteq))

raw_reaction_enthalpy = np.dot(raw_formation_enthalpies, esteq)
# print("Raw reaction enthalpy: ", raw_reaction_enthalpy)


# COPIAR LO DE ABAJO ANTES DE HACER PULL DE LO DE SALVA
# Test de la funcion reaction_enthalpies

mix5 = rd.mix.IdealGas(list_of_components)
stoichiometry = np.array([-1, -1, 1, 1])
kinetic5 = rd.Kinetics(
    mix5,
    list_of_reactions=[reaction_rate],
    stoichiometry=stoichiometry,
    kinetic_argument="concentration",
    # reaction_enthalpies=np.array([5]),
)
T = 500
P = 101325

method_reaction_enthalpy = kinetic5.reaction_enthalpies(T, P)
print(method_reaction_enthalpy)
