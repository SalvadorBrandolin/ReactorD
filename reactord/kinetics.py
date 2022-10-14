import numpy as np
from thermo.eos import R
from Mix import Abstract_Mix, IdealGas_Mix, Liquid_Mix
from Substance import Substance
from scipy.integrate import quad
from utils import vectorize


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
    stoichiometry: ndarray or list
        array or list containing the stoichiometric coefficients of
        all the substances involved in the reactive system. the
        substances that do not participate in a reaction but are present
        in the reactive system, must have a zero as stoichiometric 
        coefficient.

        Example:
        Consider the following reactions with 4 components and 1 
        inert with the substances order [A, B, C, D, I]:
        A + B --->  C + 2 D
        B + C + I ---> D + I
        the matrix of coefficients must be introduced as:

        stoichiometry = np.array[
            [-1, -1, 1, 2, 0],
            [0, -1, -1, 1, 0]
        ]

        or

        stoichiometry = [
            [-1, -1, 1, 2, 0],
            [0, -1, -1, 1, 0]
        ]

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

    def __init__(
            self, 
            mix: Abstract_Mix,
            list_of_reactions: list,
            stoichiometry: list, 
            kinetic_argument: str = 'concentration',
            **options
        ) -> None:

        
        self.reactions = list_of_reactions
        self.mix = mix
        self.kinetic_argument = kinetic_argument.lower()
        
        # Data validation               
        if np.ndim(stoichiometry) == 1:
            self.num_reactions = 1
            self.num_substances = np.shape(stoichiometry)[0]
        else:
            self.num_reactions, self.num_substances = np.shape(stoichiometry) 

        if self.num_reactions != len(list_of_reactions):
            raise IndexError(
                "'stoichiometry' rows number must be equal to" 
                " list_of_reactions' length" 
            )

        if len(mix) != self.num_substances:
            raise IndexError(
                "'stoichiometry' columns number must be equal to substances" 
                "number in 'mix' object" 
            )
        
        if self.kinetic_argument == 'concentration':
            self._composition_calculator = self.mix.concentrations
        elif self.kinetic_argument == 'partial_pressure':
            self._composition_calculator = self.mix.partial_pressures
        else:
            raise ValueError(
                f"{self.kinetic_argument} is not a valid kinetic argument"
            )

        self.stoichiometry = np.array(stoichiometry).reshape(
            self.num_reactions, self.num_substances
        )
        self.std_reaction_enthalpies = self.std_reaction_enthalpies_calc()


    @vectorize(signature='(n),(),()->(n),(m)', excluded={0})
    def kinetic_eval(
        self, 
        moles: list, 
        temperature: float, 
        pressure: float
    ) -> np.ndarray:
        """Method that evaluates the reaction rate for the reaction and 
        for the mix components. 

        Parameters
        ----------
        moles : ndarray or list
            moles of each substance
        temperature : float
            Temperature [K]
        pressure : float
            Pressure [Pa]

        Returns
        -------
        ndarray, ndarray
            
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

    def std_reaction_enthalpies_calc(self) -> np.ndarray:
        """Evaluation of the standard reaction enthalpies with the given 
        stoichiometry matrix and the formation enthalphies of the
        substances in mixture.

        Returns
        -------
        ndarray[float]
            Reaction enthalpies [J/mol]
        """
        return np.dot(self.stoichiometry, self.mix.formation_enthalpies)

""" def reaction_enthalpies(self, temperature, pressure):
        
        corr_enthalpies = 
        dh = np.dot(self.stoichiometry, cp_dt_integrals)
        
        return (dh + self.std_reaction_enthalpies) """


""" # TEST FOR A SYSTEM COMPRISED OF 2 DIFFERENT REACTIONS WITH
# 3 COMPONENTS
# A -> B reaction1
# A -> C reaction2
# [[-1,1,0],[-1,0,1]] La matriz que ingreso por teclado para estas reacciones
def reaction1(concentrations, temperature):
    return 10 * np.exp(-5000 / (R*temperature)) * concentrations[0]

def reaction2(concentrations, temperature):
    return 8 * np.exp(-6000 / (R*temperature)) * concentrations[0]

list_of_reactions = [reaction1, reaction2]

water = Substance.from_thermo_database("water")
ether = Substance.from_thermo_database("ether")
methane = Substance.from_thermo_database("methane")
mix_p = Liquid_Mix([water, ether, methane])
stoichiometry_p = ([[-1,1,0],[-1,0,1]])

cinetica = Kinetics(mix_p, list_of_reactions, stoichiometry_p)

rates_i, rate_rxns = cinetica.kinetic_eval(np.array([1, 1, 2]), 300, 101325)

print(f"velocidades por componente: {rates_i}")
print(f"velocidades por reaccion: {rate_rxns}")

#print(f"Las entalpias de reaccion son: " 
#        f"{cinetica.reaction_enthalpies(500, 101325)}")

print("\nRevision de la matriz estequiometrica")
print(f"Numero de reacciones: {cinetica.num_reactions}\n"
      f"Numero de componentes: {cinetica.total_substances}")


# ENTHALPY TEST EVALUATION FOR AN EXOTHERMIC REACTION:
# CH4 + 2 O2 ----> 2 H2O + CO2
def combustion(concentrations, temperature):
    return 10 * np.exp(-5000 / (R*temperature)) * concentrations[0]

reaction_list= [combustion]

# Substances involved:
methane = Substance.from_thermo_database("methane")
oxygen = Substance.from_thermo_database("oxygen")
water = Substance.from_thermo_database("water")
co2 = Substance.from_thermo_database("co2")

# Objects from classes Mix, Stoichiometry and Kinetics are created
exothermic_reaction_mix = IdealGas_Mix([methane, oxygen, water, co2])
stoichiometry_exothermic = ([-1, -2, 2, 1])
 
exothermic_reaction_kinetics = Kinetics(
    exothermic_reaction_mix,reaction_list,stoichiometry_exothermic)

#enthalpy = exothermic_reaction_kinetics.reaction_enthalpies(298, 101325)

#if float(enthalpy) <= 0:
#    print("\nEnthalpy is a negative value, as expected")
#else:
#    print("\nEnthalpy is not a negative number!!")

#print(f"The enthalpy of the combustion reaction is: "
#      f"{round(float(enthalpy),2)} J/mol")

#print("\nRevision de la matriz estequiometrica")
#print(f"Numero de reacciones: {exothermic_reaction_kinetics.num_reactions}\n"
#      f"Numero de componentes: {exothermic_reaction_kinetics.total_substances}") """