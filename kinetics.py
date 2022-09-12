import numpy as np

class Kinetics:
    """Kinetics object class"

    Parameters
    ----------
    ri: array of all kinetics reactions laws in the system.
    stoichiometry: array of coefficients of reagents with shape (r,c)
    Where r is the number of reactions and c number of compounds
    Each reaction row must start with the reaget coefficient which government the kinetics law.
    """


    def __init__(self, ri, stoichiometry):
        self.ri=ri
        self.stoichiometry= stoichiometry

        self.r_compounds= np.zeros(self.stoichiometry.shape)

        for i in self.stoichiometry:
            self.r_compounds[i]=self.stoichiometry[i]*self.ri(i)

#%%
ri=np.array([2,2])
st=np.array([[-1, 1, 0],[-1, 2, 3]])

a=Kinetics(ri,st)

print (a)
