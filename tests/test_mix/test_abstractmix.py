"""import numpy as np

import reactord as rd

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance

# NO SUPE COMO HACER LAS PRUEBAS DE NOT IMPLEMENTED ERROR
# PARA LA CLASE ABSTRACTA AbstractMix. LO QUE ESTÁ ABAJO ES UNA
# PRUEBA BASADA EN LO QUE PUDE VER EN:
# https://clamytoe.github.io/articles/2020/Mar/12/testing-abcs-with-abstract-methods-with-pytest/


def test_abstractmix():
    AbstractMix.__abstractmethods__ = set()

    @dataclass
    class Dummy(AbstractMix):
        pass

    test_mix = test_abstractmix()
    d = Dummy(test_mix)

    concentrations = d.concentrations(moles, temperature, pressure)
    volume = d.volume()
    formation_enthalpies_correction = d.formation_enthalpies_correction(
        temperature, pressure
    )
    formation_enthalpies_set = d._formation_enthalpies_set()
    mix_heat_capacity = d.mix_heat_capacity(moles, temperature, pressure)
"""
