"""Module to construct stoichiometry matrix."""
from typing import List

import numpy as np
from numpy.typing import NDArray

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic

from sympy import linear_eq_to_matrix, symbols


def stoichiometry_matrix_builder(
    mix: AbstractMix, equations: List[Symbolic]
) -> NDArray[np.float64]:
    """Build stoichiometry matrix.

    Parameters
    ----------
    mix : AbstractMix
    equations : List[Symbolic]


    Returns
    -------
    ndarray
        Stoichiometry matrix
    """
    n_reactions = np.size(equations)
    n_substances = len(mix)

    matrix = np.zeros((n_reactions, n_substances))

    for i, eq in enumerate(equations):
        syms = [symbols(name) for name in eq._names]
        coeffs = (
            np.array(linear_eq_to_matrix(eq._expression, syms)[0])
            .astype(np.float64)
            .ravel()
        )

        for k, name in enumerate(eq._names):
            matrix_index = np.argwhere(name == mix.names).ravel()[0]
            matrix[i, matrix_index] = coeffs[k]

    return matrix
