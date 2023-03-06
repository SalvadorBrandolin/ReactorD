import numpy as np
from numpy.typing import NDArray

def isobaric_pressure_balance(pfr) -> NDArray:
        return np.zeros(np.size(pfr.length_coordinate))