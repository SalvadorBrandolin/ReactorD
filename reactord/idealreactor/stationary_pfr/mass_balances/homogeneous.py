import numpy as np
from numpy.typing import NDArray

def homogeneous_mass_balance(pfr) -> NDArray:
        return np.multiply(pfr.substances_reaction_rates, pfr.transversal_area)
