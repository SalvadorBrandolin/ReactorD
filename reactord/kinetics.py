import numpy as np
from thermo.eos import R
from Mix import Mix
from Stoichiometry import Stoichiometry
from Substance import Substance
from scipy.integrate import quad


class Kinetics:
    """Kinetic object builder

    Parameters
    ----------
    list_of_reactions : ndarray or list [function]
        array that constains user defined python functions with the 
        form: function(composition, temperature) where 
        composition is a (number_of_components) dimension array
        that contains the partial pressures [Pa] or the concentrations
        of the substances.            
    mix : Mix object
        Mix object defined with all the substances present in the system
    stoichiometry : Stoichiometry object
        Stoichiometry object
    kinetic_argument : string
        string that indicates on wich concentration unit meassure the
        kinetic rate function are evaluated. Avaliable options: 
        'concentration', 'partial_pressure'
    enthalpy_of_reaction : ndarray or list, optional
        array that contains the enthalpy of reaction of each reaction 
        in list_of_reactions [j/mol/K]. Elements of the list may be set
        as None, then, that values are calculated using the heat
        capacity and enthalpy of formation of the substances on mix. 
        Single None value is accepted and all the values are 
        calculated, default = None.
    """

    def __init__(self, 
            list_of_reactions, 
            mix, 
            stoichiometry, 
            kinetic_argument = 'concentration',
            enthalpy_of_reaction = None
        ):
                                
        self.reactions = list_of_reactions
        self.mix = mix
        self.stoichiometry = stoichiometry.coefficients
        self.argument = kinetic_argument.lower()
        self.standard_reactions_enthalpy = self._standard_reaction_enthalpy()
        
        # The method <concentrations> from the class Mix is
        # assigned to self._composition_calculator 
        if self.argument == 'concentration':
            self._composition_calculator = self.mix.concentrations
        elif self.argument == 'partial_pressure':
            self._composition_calculator = self.mix.partial_pressures
        else:
            raise ValueError(
                f"{self.argument} is not a valid kinetic argument"
            )
        
        if enthalpy_of_reaction is None:
            self._enthalpy_of_reaction_data = np.full(
                len(self.stoichiometry), None
            )
        else:
            self._enthalpy_of_reaction_data = enthalpy_of_reaction

    def kinetic_eval(self, moles, temperature, pressure):
        """Method that evaluates the reaction 

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
            Total Pressure [Pa]

        Returns
        -------
        ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡NO SE DOCUMENTAR ESTA SALIDA!!!!!!!!!!!!!!!!!

        """
        # The partial pressures or concentrations are calculated:
        composition = self._composition_calculator(moles, temperature, pressure)
        
        # Rates for each individual reaction:
        reaction_rates = np.array(
            [reaction(composition, temperature) for reaction in self.reactions]               
        )
        # Rates for each compound:     
        rates_i = np.matmul(reaction_rates, self.stoichiometry)
        return rates_i, reaction_rates 

    def _standard_reaction_enthalpy(self):
        return np.dot(self.stoichiometry, self.mix.h_formations)

    def reaction_enthalpy(self, temperature, pressure):
        t_0 = 298.15

        if self.mix.phase == 'liquid':
            cp_dt_integrals = np.array([])

            for substance in self.mix.substances:
                cp_dt_integral, error = quad(
                    substance.heat_capacity_liquid, t_0, temperature
                )
                cp_dt_integrals = np.append(cp_dt_integrals, cp_dt_integral)
            
            dh = np.dot(self.stoichiometry, cp_dt_integrals)

        if self.mix.phase == 'gas':
            cp_dt_integrals = np.array([])

            for substance in self.mix.substances:
                cp_dt_integral, error = quad(
                    substance.heat_capacity_gas, t_0, temperature
                )
                cp_dt_integrals = np.append(cp_dt_integrals, cp_dt_integral)
            
            dh = np.dot(self.stoichiometry, cp_dt_integrals)

        return (dh + self.standard_reactions_enthalpy)

"""
#DE ACA PARA ABAJO SOLAMENTE SE PRUEBA LA CLASE:
#Ejemplo
#A -> B reaction1
#A -> C reaction2
# [[-1,1,0],[-1,0,1]] La matiz que ingreso por teclado para estas reacciones
def reaction1(concentrations, temperature):
    return 10 * np.exp(-5000 / (R*temperature)) * concentrations[0]

def reaction2(concentrations, temperature):
    return 8 * np.exp(-6000 / (R*temperature)) * concentrations[0]

list_of_reactions = [reaction1, reaction2]

water = Substance.from_thermo_database("water")
ether = Substance.from_thermo_database("ether")
methane = Substance.from_thermo_database("methane")

mix_p = Mix([water,ether, methane], "liquid")
stoichiometry_p = Stoichiometry([[-1,1,0],[-1,0,1]])

cinetica = Kinetics(list_of_reactions, mix_p, stoichiometry_p)

a, b = cinetica.kinetic_eval([1, 1, 2], 300, 101325)

print(a)
print(b)

print (f"Las entalpias de reaccion son: {cinetica.reaction_enthalpy()}")

## ALGUNAS PRUEBAS PARA EL METODO reaction_enthalpy
#A=cinetica.stoichiometry
#dot1= np.dot(A[0],cinetica.mix.formation_enthalpies)
#dot2= np.dot(A[1],cinetica.mix.formation_enthalpies)
#print(A[0], cinetica.mix.formation_enthalpies, dot1)
#print(A[1], cinetica.mix.formation_enthalpies, dot2)
"""