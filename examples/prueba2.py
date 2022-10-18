import matplotlib.pyplot as plt
import numpy as np

import reactord as rd

water = rd.Substance.from_thermo_database("water")
ethanol = rd.Substance.from_thermo_database("ethanol")

mix = rd.mix.IdealSolution([water, ethanol])

# stoichiometry definition
stoichiometry = np.array([-1, 1])

# kinetic law


def christ_reaction(concentrations, t):
    ra = 100 * np.exp(-30000 / (8.314 * t)) * concentrations[0]
    return ra


# Reactor
r_dim = np.array([0, 0.5])
area = 1
f_in = np.array([10, 0])
f_out = np.array([np.nan, np.nan])

pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
    mix=mix,
    list_of_reactions=[christ_reaction],
    stoichiometry=stoichiometry,
    reactor_dims_minmax=r_dim,
    transversal_area=area,
    pressure=101325,
    reactor_isothermic_temperature=278.15,
    reactor_f_in=f_in,
    reactor_f_out=f_out,
    kinetic_argument="concentration",
)

solution = pfr.simulate(1000)

x = solution.x

Fa, Fb = solution.y

plt.figure(0)
plt.plot(x, Fa)
plt.plot(x, Fb)

# %%
