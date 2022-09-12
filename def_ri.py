import numpy as np
from thermo.eos import R
from Mix import Mix
from Stoichiometry import Stoichiometry
from Substance import Substance

#A -> B reaction1
#A -> C reaction2

def reaction1(concentrations,T):
    return 10*np.exp(-5000/(R*T))*concentrations[0]

def reaction2(concentrations,T):
    return -8*np.exp(-6000/(R*T))*concentrations[0]


list_of_r = [reaction1, reaction2]
water= Substance.from_thermo_database("water")
ether= Substance.from_thermo_database("ether")
methane= Substance.from_thermo_database("methane")

mix_p=Mix([water,ether, methane], "liquid")
stoichiometry_p=Stoichiometry(2, 3)

class Kinetics():
    def __init__(self,list_of_r, Mix, Stoichiometry):
        self.list_of_r = list_of_r
        self.mix = Mix
        self.stoichiometry = Stoichiometry

    def kinetic_eval(self, moles,T,P):
        concentrations = self.mix.concentrations(moles,T,P)
        kinetic_evaluated = np.array([kinetic(concentrations,T) for kinetic in
                                     self.list_of_r])
        kinetic_stoichiometry = np.array([ self.stoichiometry[reaction]*kinetic_evaluated[reaction] 
                                         for reaction in self.stoichiometry])
        return kinetic_stoichiometry

cinetica = Kinetics(list_of_r,mix_p,stoichiometry_p)

print (cinetica.kinetic_eval([1, 1, 2], 300, 1))