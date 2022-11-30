import numpy as np
import matplotlib.pyplot as plt
# import the reactord package
import reactord as rd

# The reaction is as follows:
# A ---> B

# The substances objects are instantiated:
A = rd.Substance.from_thermo_database("cis-2-butene")
B = rd.Substance.from_thermo_database("trans-2-butene")

# The names of the substances can be changed for convenience:
A.name = "A"
B.name = "B"

# The mixture object is created:
mixture = rd.mix.IdealSolution([A, B])

# A function for the reaction rate is defined:
def rate(concentrations, temperature):
    return 0.00001*concentrations[0]



fa_in = 5 / 3600  # mol/s

f_volumetric = 10 * 0.001 / 3600  # m3/s

v_pfr = 99 * 0.001  # m3


# The tubular reactor is instantiated:
pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
    mix=mixture,
    list_of_reactions=[rate],
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

largo_reactor = solution.x
var_dependientes = solution.y

fig1, ax1 = plt.subplots(1)
ax1.plot(largo_reactor, var_dependientes[0], "-r", label="Flux_A", linewidth=3)
ax1.plot(largo_reactor, var_dependientes[1], "-g", label="Flux_B", linewidth=3)
#ax1.set(xlabel="Reactor volume (m^3)", ylabel="molar fluxes (mol/s)")
ax1.set_xlabel("Reactor volume (m^3)", fontsize=10)
ax1.set_ylabel("molar fluxes (mol/s)", fontsize=10)

ax1.legend()
fig1.savefig("Butene_Isomerization.png")

plt.show()
