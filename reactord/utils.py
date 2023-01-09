"""Utils module.

Utility code for ReactorD library.
"""
from typing import Callable

import numpy as np


def vectorize(signature: str = None, excluded: set = None) -> Callable:
    """Vectorize user_function.

    Decorator that applies the numpy.vectorize to a given user_function.

    Parameters
    ----------
    signature : str, optional
        Input and output dimensions, by default None
        Example:
        '(m,n),(n)->(m)' for vectorized matrix-vector multiplication
    excluded : set, optional
        Set of strings or integers representing the positional or
        keyword arguments for which the function won't be vectorized.
        These will be passed directly to user_function unmodified.
        By default None

    Returns
    -------
    Callable
        Array with given dimensions for evaluated user_function
    """

    def outer(user_function: Callable) -> Callable:
        def inner(*args, **kwargs):
            f = np.vectorize(
                user_function, signature=signature, excluded=excluded
            )
            return f(*args, **kwargs)

        return inner

    return outer
