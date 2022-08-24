import numpy as np

class Stoichiometry:
    def __init__(self, reactions, total_comp, names=None):
        self.reactions=reactions
        self.total_comp= total_comp
        self.names= names #A label with the reaction name

        #The stoichiometry array is inicialized:
        self.matrix = np.zeros((reactions, total_comp))

        self.matrix = input (f"Introduce the matrix of coefficients with"
f" {self.reactions} files and {self.total_comp} columns."
f"\nUse negative numbers for reactants and positive for products:")
        
   
    def __str__ (self):
        return ((f"Total components: {str(self.total_comp)}\n") + 
        (f"Number of reactions: {str(self.reactions)} \n"))


# ACA SE PRUEBA EL OBJETO    
prueba= Stoichiometry(2,2,["reaccion1", "esterficacion"])
print (prueba.matrix) #La matriz de coefficientes
print (prueba.names) #Nombres de las reacciones
print(prueba) #Se prueba el metodo __str__
