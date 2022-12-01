import numpy as np

import reactord as rd

"""Fogler fourth ed. P1.15a as an initial value problem"""


def volume(temperature, pressure):
    return 1 / (fa_in / f_volumetric)


def kinetic(concentrations, temperature):
    return k


# Fogler's exact solution


def fogler(vol):
    concentration = (fa_in / f_volumetric) - vol * k / f_volumetric

    return concentration


fogler_concentrations = np.array([])

fa_in = 5 / 3600  # mol/s

f_volumetric = 10 * 0.001 / 60  # m3/s

k = 0.05 / 3600 / 0.001  # mol/s/m3

v_pfr = 99 * 0.001  # m3

substance_a = rd.Substance(name="A")
substance_a.volume_liquid = volume

substance_b = rd.Substance(name="B")
substance_b.volume_liquid = volume

mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
    mix=mixture,
    list_of_reactions=[kinetic],
    stoichiometry=[-1, 1],
    kinetic_argument="concentration",
    reactor_dim_minmax=[0, v_pfr],
    transversal_area=1,
    isothermic_temperature=298.15,
    isobaric_pressure=101325,
    molar_flow_in={"A": fa_in, "B": 0},
)

# reactord solution
solution = pfr.simulate(grid_size=100)
reactord_concentrations = np.array([])

# Comparisson

for i, v in enumerate(solution.x):
    reactord_concentrations = np.append(
        reactord_concentrations,
        mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
    )

    fogler_concentrations = np.append(fogler_concentrations, fogler(v))

assert np.allclose(reactord_concentrations, fogler_concentrations)
