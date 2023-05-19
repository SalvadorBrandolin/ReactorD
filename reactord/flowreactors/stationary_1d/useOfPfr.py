import reactord as rd
import numpy as np
import reactord.flowreactors.stationary_1d.pfr as pfr
import matplotlib.pyplot as plt


A = rd.Substance.from_thermo_database("MeOH", "methanol")
B = rd.Substance.from_thermo_database("H2O", "H2O")

mixture = rd.mix.IdealSolution([A, B])

def r_rate(c, t, cons):
    return np.full(np.size(t), cons["k"])

reactions = {"r1": {"eq": A > B, "rate": r_rate}}
constantes_cineticas = {"k": 1e-5}

cinetica = rd.Kinetic(mixture, reactions, constantes_cineticas)

mb = pfr.mass_balances.MolarFlow(molar_flows_in={"MeOH": 100, "H2O": 0})
eb = pfr.energy_balances.Isothermic(350)
pb = pfr.pressure_balances.Ergun(
    pressure={"in": 10 * 101325}, porosity=0.45, particle_diameter=0.00635
)


reactor = pfr.PFR(
    kinetic=cinetica,
    reactor_length=10,
    transversal_area=0.2,
    grid_size=20,
    mass_balance=mb,
    energy_balance=eb,
    pressure_balance=pb,
)

results = reactor.simulate()
X = reactor.ode_solution

# print(len(reactor.ode_solution))
print(X.y)
# self.sim_df = pd.DataFrame(result.T, columns=columns, index=None)



plt.plot(X.x, X.y[0], "-g",
         X.x, X.y[1], "-b",
         X.x, X.y[2], "--r",
                 
         )




plt.show()


#x = np.linspace(0, 20, 100)
#plt.plot(x, np.sin(x))
#plt.show()