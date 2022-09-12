import numpy as np
from thermo.eos import R
from Mix import Mix
from Stoichiometry import Stoichiometry
from Substance import Substance


class Kinetics():
    """The class Kinetics returns the rate of reaction of each 
    species present in the reactor"""
    def __init__(self, list_of_r, Mix, Stoichiometry):
        self.rxns = list_of_r
        self.mix = Mix
        self.stoichiometry = Stoichiometry.coefficients

#ABAJO EXPLICO EL PORQUE DE LA FUNCION KINETIC_EVAL
    def kinetic_eval(self, moles, temp, presure, compound_rate = True):
        concentrations = self.mix.concentrations(moles, temp, presure)
        rxn_rates = np.array(
            [rxn(concentrations, temp)for rxn in self.rxns]               
        )
             
        if compound_rate == True:
            rates_i = np.matmul(rxn_rates, self.stoichiometry)
            return rates_i 
        else:
            return rxn_rates

""" A la funcion kinetic_eval, le agregué el parametro compound_rate=True.
Lo hice así porque primero quise hacer que devuelva los arrays rxn_rates y
rates_i en una sola sentencia. Esto me tiró un warning que me
dice que se desaconseja hacerlo de esa manera. Asi que para evitar eso,
escribí la funcion para que devuelva una u otra cosa. Como generalmente
nos interesa las velocidades por componente, se devuelven esas por defecto.
Si quieren ver como era antes, comenten las lineas 17 a 27 y descomenten las
39 a 49 """


"""def kinetic_eval(self, moles, temp, presure):
        concentrations = self.mix.concentrations(moles, temp, presure)
        rxn_rates = np.array(
            [rxn(concentrations, temp)for rxn in self.rxns]
        ) # rxn(concentrations,T) es cada funcion de la lista
          # de reacciones. Es decir es la velocidad bajo condiciones
          # especificas (para cada reaccion) 

        rates_i = np.matmul(rxn_rates, self.stoichiometry)

        return rates_i, rxn_rates""" # CREAR UN NDARRAY DE ESTAS
                                                # CARACTERISTICAS ESTA
                                                # DEPRECADO
 


#DE ACA PARA ABAJO SOLAMENTE SE PRUEBA LA CLASE:
#Ejemplo
#A -> B reaction1
#A -> C reaction2
# [[-1,1,0],[-1,0,1]] La matiz que ingreso por teclado para estas reacciones
def reaction1(concentrations, temp):
    return 10 * np.exp(-5000 / (R*temp)) * concentrations[0]

def reaction2(concentrations, temp):
    return 8 * np.exp(-6000 / (R*temp)) * concentrations[0]

list_of_r = [reaction1, reaction2]

water = Substance.from_thermo_database("water")
ether = Substance.from_thermo_database("ether")
methane = Substance.from_thermo_database("methane")

mix_p = Mix([water,ether, methane], "liquid")
stoichiometry_p = Stoichiometry(2, 3)

cinetica = Kinetics(list_of_r, mix_p, stoichiometry_p)

a = cinetica.kinetic_eval([1, 1, 2], 300, 101325, compound_rate=False)
print(a, np.shape(a),type(a))

b = cinetica.kinetic_eval([1, 1, 2], 300, 101325)
print(b, np.shape(b),type(b))

#PRUEBAS DE NUMPY VARIAS
#B= np.array([[-1,1,0], [-1,0,1]])
#print('\n', B, np.shape(B))
#print("fila o columna?:", B[0])


