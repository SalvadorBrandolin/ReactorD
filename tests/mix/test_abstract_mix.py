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
        NewMixture(
            substance_list=[],
            phase_nature="liquid",
            viscosity_mixing_rule="linear",
        )


def test_abastract_class_not_implemented_errors():
    class NewMixture(rd.mix.AbstractMix):
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

    mixture = NewMixture(
        substance_list=[],
        phase_nature="liquid",
        viscosity_mixing_rule="linear",
    )

    with pytest.raises(NotImplementedError):
        mixture.volume([1, 1], 298.15, 101325)

    with pytest.raises(NotImplementedError):
        mixture.mix_heat_capacity([1, 1], 298.15, 101325)

    with pytest.raises(NotImplementedError):
        mixture._formation_enthalpies_set()

    with pytest.raises(NotImplementedError):
        mixture.formation_enthalpies_correction(298.15, 101325)


def test_abastract_class_viscosity_mixing_rule():
    class AnotherMixture(rd.mix.AbstractMix):
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

    mixture = AnotherMixture(
        substance_list=[],
        phase_nature="liquid",
        viscosity_mixing_rule="linear",
    )

    assert mixture.phase_nature == "liquid"
    assert mixture.viscosity_mixing_rule == "linear"

    mixture = AnotherMixture(
        substance_list=[],
        phase_nature="gas",
        viscosity_mixing_rule="herning_zipperer",
    )

    assert mixture.phase_nature == "gas"
    assert mixture.viscosity_mixing_rule == "herning_zipperer"

    with pytest.raises(ValueError):
        AnotherMixture(
            substance_list=[],
            phase_nature="hierophant green",
            viscosity_mixing_rule="linear",
        )

    with pytest.raises(ValueError):
        AnotherMixture(
            substance_list=[],
            phase_nature="gas",
            viscosity_mixing_rule="Liu Kang",
        )
