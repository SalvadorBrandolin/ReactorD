import numpy as np
import matplotlib.pyplot as plt
import reactord as rd

# Primero describimos una reacción, por ejemplo
# la bromacion del benceno


# C6H6     +     Br2 ----> C6H5Br      +      HBr
# benzene      Bromine   bromobenzene   hydrobromic acid
#    A     +      B  ---->    C        +       D

# Supongamos ademas que tenemos una reaccion secundaria de neutralizacion:
# NaOH    +    HBr   ---->  NaBr    +     H2O
#   E     +     D    ---->   F      +      G

# Esto define la siguiente estequiometria:
stoichiometry = [[-1, -1, 1, 1, 0, 0, 0], [0, 0, 0, -1, -1, 1, 1]]

# Ahora creamos los objetos Substance que participan en ambas reacciones:

ben = rd.Substance.from_thermo_database("benzene")
br2 = rd.Substance.from_thermo_database("bromine")
brben = rd.Substance.from_thermo_database("bromobenzene")
hbr = rd.Substance.from_thermo_database("10035-10-6")  # using CAS number
naoh = rd.Substance.from_thermo_database("naoh")
nabr = rd.Substance.from_thermo_database("nabr")
h2o = rd.Substance.from_thermo_database("water")

ben.name = "A"
br2.name = "B"
brben.name = "C"
hbr.name = "D"
naoh.name = "E"
nabr.name = "F"
h2o.name = "G"

# print (water.name, water.molecular_weight)


# Comprobamos que las sustancias creadas sean los compuestos deseados
"""
print(ben.name, ben.molecular_weight, ben.normal_boiling_point)
print(br2.name, br2.molecular_weight, br2.normal_boiling_point)
print(brben.name, brben.molecular_weight, brben.normal_boiling_point)
print(hbr.name, hbr.molecular_weight, hbr.normal_boiling_point)
"""

# Ahora creamos una mezcla líquida, que contiene todas las sustancias
# involucradas. Usamos para esto la clase IdealSolution y asumimos que todo
# es una solucion ideal (aunque en realidad no es asi):
mezcla = rd.mix.IdealSolution(
    ben=ben, br2=br2, brben=brben, hbr=hbr, naoh=naoh, nabr=nabr, h2o=h2o
)

# Definimos funciones para cada una de las velocidades de reaccion:


def rate1(concentrations, temperature):
    return 0.000003 * concentrations[0] * concentrations[1]


def rate2(concentrations, temperature):
    return 0.000005 * concentrations[3] * concentrations[4]


# Definimos los flujos molares de entrada:
flu_entrada = np.asarray([10, 10, 0, 0, 2, 0, 0]) * 1000 / 3600  # moles/s
vol_pfr = 1000 * 0.0001  # m³

# print(flu_entrada, np.shape(flu_entrada))
# for i in range(np.size(flu_entrada)):
#    print(flu_entrada[i])

# Instanciamos el reactor de flujo en piston
"""pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
    mezcla, 
    list_of_reactions=[rate1, rate2],
    stoichiometry=stoichiometry,
    kinetic_argument="concentration",
    reactor_dim_minmax=[0,vol_pfr],
    transversal_area=1,
    isothermic_temperature=350, # K
    isobaric_pressure=200000, # Pa
    molar_flow_in={"ben":flu_entrada[0], "br2":flu_entrada[1],
    "brben":flu_entrada[2], "hbr":flu_entrada[3], "naoh":flu_entrada[4],
    "nabr":flu_entrada[5], "h2o":flu_entrada[6]},
    #molar_flow_out=None,
    #catalyst_particle=None,
)
"""


pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
    mezcla,
    list_of_reactions=[rate1, rate2],
    stoichiometry=stoichiometry,
    kinetic_argument="concentration",
    reactor_dim_minmax=[0, vol_pfr],
    transversal_area=1,
    isothermic_temperature=298.15,  # K
    isobaric_pressure=101325,  # Pa
    molar_flow_in={
        "A": flu_entrada[0],
        "B": flu_entrada[1],
        "C": flu_entrada[2],
        "D": flu_entrada[3],
        "E": flu_entrada[4],
        "F": flu_entrada[5],
        "G": flu_entrada[6],
    },
    # molar_flow_out=None,
    # catalyst_particle=None,
)

# print (pfr.molar_flow_in)

# Simulamos el funcionamiento del reactor
solution = pfr.simulate(grid_size=100)
reactor_volume = solution.x
pfr_concentrations = solution.y


fig, ax = plt.subplots(1)
ax.plot(reactor_volume, pfr_concentrations[0], "-r", linewidth=2, label="")
ax.plot(reactor_volume, pfr_concentrations[1], "--r", linewidth=2)
ax.plot(reactor_volume, pfr_concentrations[2], "-g", linewidth=2)
ax.plot(reactor_volume, pfr_concentrations[3], "--g", linewidth=2)

ax.set_xlabel(r"Reactor volume [$ m{³}$]")
ax.set_ylabel(r"Flujos molares [$mol/s$]")

ax.legend()

plt.show()
