from typing import Callable

import numpy as np


def vectorize(signature: str = None, excluded: set = None) -> Callable:
    def outer(user_function: Callable) -> Callable:
        def inner(*args, **kwargs):
            f = np.vectorize(
                user_function, signature=signature, excluded=excluded
            )
            return f(*args, **kwargs)

        return inner

    return outer
