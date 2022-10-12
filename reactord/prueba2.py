#%%
from ..reactord import Substance
import numpy as np

water = rd.Substance.from_thermo_database('water')
ethanol = rd.Substance.from_thermo_database('ethanol')

#%%
#mixture definition
mix = rd.Liquid_Mix([water, ethanol])

#stoichiometry definition
stoichiometry = np.array([-1,1])
#kinetic law
def christ_reaction(concentrations, T):
    ra = 100 * np.exp(-30000 / (8.314 * T)) * concentrations[0]
    return ra

#kinetic
kinetic = Kinetics(np.array([christ_reaction]), mix, stoichiometry)

#Reactor
r_dim = np.array([0, 0.5])
area = 1
f_in = np.array([10, 0])
f_out = np.array(['var', 'var'])
t_in = 298.15
t_out = 'var'

pfr = Homogeneous_PFR(mix, kinetic, r_dim, area, 101325, 'non-isothermal',
                      f_in, f_out, t_in, t_out)

#%%
solution = pfr.solve(100)

#%%
import matplotlib.pyplot as plt
x = solution.x

Fa, Fb, T, P, Ta = solution.y

plt.figure(0)
plt.plot(x, Fa)
plt.plot(x, Fb)

plt.figure(1)

plt.plot(x, T)
plt.plot(x, Ta)

plt.figure(2)
plt.plot(x, P)