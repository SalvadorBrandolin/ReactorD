import reactord as rd
import reactord.flowreactors.stationary_1d.pfr as pfr

import numpy as np

a = rd.Substance.from_thermo_database(name="methane", thermo_identification="methane")
b = rd.Substance.from_thermo_database(name="ethane", thermo_identification="ethane")

mix = rd.mix.IdealGas([a, b])

def ra(concentrations, temperature):
    ra = 0.05 / 3600 / 0.001 * concentrations[0]
    return ra

stoichiometry = np.array([-1, 1])

kinetic = rd.Kinetics(mix, [ra], stoichiometry)

mb = pfr.mass_balances.df_dz(
    molar_flows_in={"methane": 0.001, "ethane": 0},
    molar_flows_out={}
)

eb = pfr.energy_balances.Isothermic(400)

pb = pfr.pressure_balances.Isobaric(101325)

reactor = pfr.PFR(
    mix = mix,
    kinetics = kinetic,
    reactor_length = 0.01,
    transversal_area = 1,
    grid_size=100,
    mass_balance=mb,
    energy_balance=eb,
    pressure_balance=pb
)

sol = reactor.simulate()