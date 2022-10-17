import numpy as np


def vectorize(signature=None, excluded=None):
    def outer(user_function):
        def inner(*args, **kwargs):
            f = np.vectorize(
                user_function, signature=signature, excluded=excluded)
            return f(*args, **kwargs)
        return inner
    return outer