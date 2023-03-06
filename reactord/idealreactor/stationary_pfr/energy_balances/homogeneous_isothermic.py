import numpy as np
from numpy.typing import NDArray

def isothermic_energy_balance(pfr) -> NDArray:
        return np.zeros(np.size(pfr.length_coordinate))