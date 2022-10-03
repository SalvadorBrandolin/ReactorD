from Substance import Substance
from Mix import Mix
from Stoichiometry import Stoichiometry
from kinetics import Kinetics
from Homogeneous_PFR import Homogeneous_PFR
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
stoichiometry = Stoichiometry(np.array([-1,1]))
#kinetic law
def christ_reaction(concentrations, T):
    ra = 0.01 * np.exp(-30000 / (8.314 * T)) * concentrations[0]
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

solution = pfr.solve(100)