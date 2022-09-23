import numpy as np


class Stoichiometry:
    """Stoichiometry object class.
    
    The stoichiometric coefficients of the reactions must be
    introduced as a numpy array. Reactants are introduced as negative
    numbers, products as positive numbers and inerts as zero.

    Example:
    Consider the following reactions with 4 components and 1 inert with 
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
        array or list containing the stoichiometric coefficients of
        all the substances involved in the reactive system. the
        substances that do not participate in a reaction but are present
        in the reactive system, must have a zero as stoichiometric 
        coefficient.
    ----------
    """
    
    def __init__(self, stoichiometric_coefficients):
        self.coefficients = stoichiometric_coefficients
        
        # num_reactions and total_comp are computed differently in
        # systems with one reaction or more than one reaction        
        if len(np.shape(stoichiometric_coefficients)) == 1:
            self.num_reactions = 1
            self.total_comp = np.shape(stoichiometric_coefficients)[0]
        else:
            self.num_reactions, self.total_comp = np.shape(
            stoichiometric_coefficients) 

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
    
""" # ACA SE PRUEBA LA CLASE    
one_reaction = Stoichiometry([-1,1,0])
two_reactions = Stoichiometry ([[1,0,2],[-1,2,3]])
print (f"Una reaccion: {one_reaction}")
print (f"Dos reacciones: {two_reactions}") """
