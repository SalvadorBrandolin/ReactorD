import numpy as np

import reactord as rd


def test_vectorization():
    @rd.vectorize(signature="(n),()->()", excluded={0})
    def vectorized_function(ignore_arg, vector, scalar):
        sumatory = np.sum(vector)
        return scalar * sumatory

    vectors = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]])

    scalars = np.array([2, 3, 4, 5])

    ignored = np.zeros(50)

    vectorized_result = vectorized_function(ignored, vectors, scalars)

    assert np.allclose(vectorized_result, np.array([6, 18, 36, 60]))
