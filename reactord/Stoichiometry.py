import numpy as np


class Stoichiometry:
    """Stoichiometry object class.
    
    The creation of an object from this class prompts the user to
    introduce the parameters. After that, the stoichiometric
    coefficients of the reactions must be introduced as a numpy
    array. Reactants are introduced as negative numbers, products
    as positive number and inerts as zero.

    Example:
    Consider the following reactions with 4 components and 1 inert and 
    the substances order [A, B, C, D, I]:
    A + B --->  C + 2 D
    B + C + I ---> D + I
    the matrix of coefficients must be introduced as:

    stoichiometric_coefficients = np.array[
        [-1, -1, 1, 2, 0],
        [0, -1, -1, 1, 0]
    ]

    or

    stoichiometric_coefficients = [
        [-1, -1, 1, 2, 0],
        [0, -1, -1, 1, 0]
    ]

    Parameters
    ----------
    stoichiometric_coefficients: ndarray or list
        array or list containing all the stoichiometric coefficients of
        all the substances involved in the reactive system. the
        substances that do not participate in a reaction but are present
        in the reactive system, it must have a zero as
        stoichiometric coefficient.
    ----------
    """
    
    def __init__(self, stoichiometric_coefficients):
        self.num_reactions, self.total_comp = np.shape(
            stoichiometric_coefficients) 
        self.coefficients = stoichiometric_coefficients
              
    def __str__(self):
        return(f"Total components: {str(self.total_comp)}\n"
        f"Number of reactions: {str(self.num_reactions)} \n")

    def __iter__(self):
      return self.coefficients.__iter__()
    
    def __len__(self):
        return len(self.coefficients)
    
    #selecctor de elemento dentro del objeto
    #def __getitem__(self, index):
    #    return self.coefficients[index]
    #buscador de elemento detro del objeto
    
# ACA SE PRUEBA LA CLASE    
#prueba= Stoichiometry(2,2)
#print(prueba.coefficients,np.shape(prueba.coefficients)) #La matriz de coeficientes
#print(prueba.names) #Nombres de las reacciones
#print(prueba) #Se prueba el metodo __str__