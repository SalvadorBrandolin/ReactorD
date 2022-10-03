import numpy as np

def vectorize(signature=None):
    def outer(user_function):
        def inner(*args, **kwargs):
            f = np.vectorize(user_function, signature=signature)
            return f(*args, **kwargs)
        return inner
    return outer


@vectorize(signature='(n)->(n)')
def func(x):
    return 2*x

a = [1,2,3,4,5]


@vectorize(signature=None)
def func2(x,y):
    return x+y

print(func2(a,a))