from typing import List

import pytest

import reactord as rd


def test_abstract_class_type_error():
    with pytest.raises(TypeError):
        rd.mix.AbstractMix()


def test_not_defining_abastract_methods():
    class NewMixture(rd.mix.AbstractMix):
        ...

    with pytest.raises(TypeError):
        mixture = NewMixture()


def test_abastract_class_not_implemented_errors():
    class NewMixture(rd.mix.AbstractMix):
        def concentrations(
            self, moles: List[float], temperature: float, pressure: float
        ):
            return super().concentrations(moles, temperature, pressure)

        def volume(
            self, moles: List[float], temperature: float, pressure: float
        ):
            return super().volume(moles, temperature, pressure)

        def mix_heat_capacity(
            self, moles: List[float], temperature: float, pressure: float
        ):
            return super().mix_heat_capacity(moles, temperature, pressure)

        def _formation_enthalpies_set(self):
            return super()._formation_enthalpies_set()

        def formation_enthalpies_correction(
            self,
            temperature: float,
            pressure: float,
        ):
            return super().formation_enthalpies_correction(
                temperature, pressure
            )

    mixture = NewMixture()

    with pytest.raises(NotImplementedError):
        mixture.concentrations([1, 1], 298.15, 101325)

    with pytest.raises(NotImplementedError):
        mixture.volume([1, 1], 298.15, 101325)

    with pytest.raises(NotImplementedError):
        mixture.mix_heat_capacity([1, 1], 298.15, 101325)

    with pytest.raises(NotImplementedError):
        mixture._formation_enthalpies_set()

    with pytest.raises(NotImplementedError):
        mixture.formation_enthalpies_correction(298.15, 101325)
