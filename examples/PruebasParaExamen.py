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
    return 3e-6 * concentrations[0] * concentrations[1]


def rate2(concentrations, temperature):
    return 1e-5 * concentrations[3] * concentrations[4]


# Definimos los flujos molares de entrada:
flu_entrada = np.asarray([10, 13, 0, 0, 6, 0, 0]) * 1000 / 3600  # moles/s

# Instanciamos el reactor de flujo en piston. Hacemos graficas para dos simula
# ciones diferentes, aumentando el volumen del reactor

volume = np.asarray([225, 230]) * 0.001
# NOTA: Al aumentar el volumen arriba de 238 *0.001, la simulacion comienza
# a tardar un huevo y las graficas salen feas.
# Parece ser que el detonante de este comportamiento es que uno de los flujos de
# reactivos se hace negativo, aunque no se por qué.
num_col = 0
plots = 2
fig, ax = plt.subplots(1, plots)

for v in volume:

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
        mezcla,
        list_of_reactions=[rate1, rate2],
        stoichiometry=stoichiometry,
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v],
        transversal_area=1,
        isothermic_temperature=350,  # K
        isobaric_pressure=200000,  # Pa
        molar_flow_in={
            ben.name: flu_entrada[0],
            br2.name: flu_entrada[1],
            brben.name: flu_entrada[2],
            hbr.name: flu_entrada[3],
            naoh.name: flu_entrada[4],
            nabr.name: flu_entrada[5],
            h2o.name: flu_entrada[6],
        },
        # molar_flow_out=None,
        # catalyst_particle=None,
    )

    # Simulamos el funcionamiento del reactor
    solution = pfr.simulate(grid_size=50)
    reactor_volume = solution.x
    pfr_concentrations = solution.y


    ax[num_col].plot(
        reactor_volume,
        pfr_concentrations[0],
        "-r",
        linewidth=2,
        label="benzene",
    )
    ax[num_col].plot(
        reactor_volume, pfr_concentrations[1], "--r", linewidth=2, label="Br2"
    )
    ax[num_col].plot(
        reactor_volume,
        pfr_concentrations[2],
        "-g",
        linewidth=2,
        label="Bromobenzene",
    )
    ax[num_col].plot(
        reactor_volume, pfr_concentrations[3], "-k", linewidth=2, label="HBr"
    )
    ax[num_col].plot(
        reactor_volume, pfr_concentrations[4], "--g", linewidth=2, label="NaOH"
    )
    ax[num_col].plot(
        reactor_volume, pfr_concentrations[5], "--b", linewidth=2, label="NaBr"
    )
    ax[num_col].plot(
        reactor_volume, pfr_concentrations[6], "-b", linewidth=2, label="H2O"
    )

    num_col += 1
    
plt.show()
