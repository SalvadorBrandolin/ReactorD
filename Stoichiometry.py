import numpy as np


class Stoichiometry:
    """Stoichiometry object class.

    Parameters
    ----------
    reactions: total number of reactions ocurring simultaneously
    total_comp: number of total components (include inerts).
    ----------

    The creation of an object from this class prompts the user to
    introduce the parameters. After that, the stoichiometric
    coefficients of the reactions must be introduced as a numpy
    array. Reactants are introduced as negative numbers, products
    as positive number and inerts as zero.

    Example:
    Consider the following reactions with 4 components and 1 inert
    A + B --->  C + 2 D
    B + C + I ---> D + I
    the matrix of coefficients must be introduced as:

    [[-1 -1 1 2 0] [0 -1 -1 1 0]]
    """
    
    def __init__(self, reactions, total_comp, names=None):
        self.reactions=reactions
        self.total_comp= total_comp
        self.names= names 

        #The coefficients array is created:        
        coefficients = input(f"Introduce the matrix of coefficients with "
            f"{self.reactions} files and {self.total_comp} columns."
            f"\nUse negative numbers for reactants and positive for products:")
        self.coefficients= np.array(coefficients)
           
    def __str__(self):
        return(f"Total components: {str(self.total_comp)}\n"
        f"Number of reactions: {str(self.reactions)} \n")

    def __iter__(self):
      return self.coefficients.__iter__()
    
    def __len__(self):
        return len(self.coefficients)
    #selecctor de elemento dentro del objeto
    def __getitem__(self, index):
        return self.coefficients[index]
    #buscador de elemento detro del objeto
    


# ACA SE PRUEBA LA CLASE    
prueba= Stoichiometry(2,2)
print(type(prueba.coefficients),np.shape(prueba.coefficients)) #La matriz de coeficientes
#print(prueba.names) #Nombres de las reacciones
#print(prueba) #Se prueba el metodo __str__
