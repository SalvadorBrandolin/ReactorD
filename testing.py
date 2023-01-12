import numpy as np

def funcion(x):
    return x+2

funcionv = np.vectorize(funcion, signature='()->()')


print(funcionv(15))