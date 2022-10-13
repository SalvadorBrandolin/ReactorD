#%%
from substance import Substance
from Mix import Liquid_Mix, IdealGas_Mix
from kinetics import Kinetics
from reactord.pfr_homogeneous_stationary import Homogeneous_PFR
import numpy as np
""" Reaccion quimica agua -> etanol, con cinetica kCagua.
Al reactor ingresa 10 mol/tiempo de agua pura a 298.15 K
La reaccions e realiza a pression atmosférica. Y es adiabatica. """

#Substance definition
water = Substance.from_thermo_database('water')
ethanol = Substance.from_thermo_database('ethanol')

#mixture definition
mix = Mix(np.array([water, ethanol]), 'liquid')

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

# %%
