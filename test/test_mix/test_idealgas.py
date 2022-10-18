import numpy as np

import reactord as rd


def test_one_substance_mix():
    hexane = rd.Substance.from_thermo_database('methane')

    oxygen = rd.Substance.from_thermo_database('oxygen')

    hydrogen = rd.Substance.from_thermo_database('hydrogen')


    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000]) 