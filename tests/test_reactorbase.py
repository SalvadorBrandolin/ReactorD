import numpy as np

import pytest

import reactord as rd


def test_abstract_error_raises():
    with pytest.raises(TypeError):
        a = rd.Substance()
        b = rd.Substance()

        mixture = rd.mix.IdealSolution(A=a, B=b)

        rd.ReactorBase(
            mix=mixture,
            list_of_reactions=[],
            stoichiometry=[-1, 1],
            kinetic_argument="concentration",
        )


def test_heritage_methods():
    # =========================================================================
    # Not implemented all abstract methods
    # =========================================================================
    class SpecificReactor(rd.ReactorBase):
        pass

    with pytest.raises(TypeError):
        SpecificReactor()


def test_not_implemented_error():
    "Test for not implemented error of reactor's intherface methods."

    class SpecificReactor(rd.ReactorBase):
        @classmethod
        def from_isothermic_isobaric(cls) -> None:
            return super().from_isothermic_isobaric()

        @classmethod
        def from_isothermic_noisobaric(cls) -> None:
            return super().from_isothermic_noisobaric()

        @classmethod
        def from_adiabatic_isobaric(cls) -> None:
            return super().from_adiabatic_isobaric()

        @classmethod
        def from_adiabatic_noisobaric(cls) -> None:
            return super().from_adiabatic_noisobaric()

        @classmethod
        def from_noisothermic_isobaric(cls) -> None:
            return super().from_noisothermic_isobaric()

        @classmethod
        def from_noisothermic_noisobaric(cls) -> None:
            return super().from_noisothermic_noisobaric()

        def _set_catalyst_operation(self):
            return super()._set_catalyst_operation()

        def _set_thermal_operation(self):
            return super()._set_thermal_operation()

        def _set_pressure_operation(self):
            return super()._set_pressure_operation()

        def _grid_builder(self) -> None:
            return super()._grid_builder()

        def _border_cond_and_initial_guesses(self) -> None:
            return super()._border_cond_and_initial_guesses()

        def _refrigerant_energy_balance(self) -> None:
            return super()._refrigerant_energy_balance()

        def simulate(self) -> None:
            return super().simulate()

        def _homogeneous_mass_balance(self) -> None:
            return super()._homogeneous_mass_balance()

        def _heterogeneous_mass_balance(self) -> None:
            return super()._heterogeneous_mass_balance()

        def _isothermic_energy_balance(self) -> None:
            return super()._isothermic_energy_balance()

        def _homogeneous_adiabatic_energy_balance(self) -> None:
            return super()._homogeneous_adiabatic_energy_balance()

        def _heterogeneous_adiabatic_energy_balance(self) -> None:
            return super()._heterogeneous_adiabatic_energy_balance()

        def _homogeneous_non_isothermic_energy_balance(self) -> None:
            return super()._homogeneous_non_isothermic_energy_balance()

        def _heterogeneous_non_isothermic_energy_balance(self) -> None:
            return super()._heterogeneous_non_isothermic_energy_balance()

        def _isobaric_pressure_balance(self) -> None:
            return super()._isobaric_pressure_balance()

        def _non_isobaric_pressure_balance_packed_bed_reactor(self) -> None:
            return super()._non_isobaric_pressure_balance_packed_bed_reactor()

        def _non_isobaric_pressure_balance_gas_phase_reaction(self) -> None:
            return super()._non_isobaric_pressure_balance_gas_phase_reaction()

        def _homogeneous_solver(self) -> None:
            return super()._homogeneous_solver()

        def _heterogeneous_solver(self) -> None:
            return super()._heterogeneous_solver()

    reactor = SpecificReactor()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_isothermic_isobaric()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_isothermic_noisobaric()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_adiabatic_isobaric()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_adiabatic_noisobaric()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_noisothermic_isobaric()

    with pytest.raises(NotImplementedError):
        SpecificReactor.from_noisothermic_noisobaric()

    with pytest.raises(NotImplementedError):
        reactor._set_catalyst_operation()

    with pytest.raises(NotImplementedError):
        reactor._set_thermal_operation()

    with pytest.raises(NotImplementedError):
        reactor._set_pressure_operation()

    with pytest.raises(NotImplementedError):
        reactor._grid_builder()

    with pytest.raises(NotImplementedError):
        reactor._border_cond_and_initial_guesses()

    with pytest.raises(NotImplementedError):
        reactor._refrigerant_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor.simulate()

    with pytest.raises(NotImplementedError):
        reactor._homogeneous_mass_balance()

    with pytest.raises(NotImplementedError):
        reactor._heterogeneous_mass_balance()

    with pytest.raises(NotImplementedError):
        reactor._isothermic_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._homogeneous_adiabatic_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._heterogeneous_adiabatic_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._homogeneous_non_isothermic_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._heterogeneous_non_isothermic_energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._isobaric_pressure_balance()

    with pytest.raises(NotImplementedError):
        reactor._non_isobaric_pressure_balance_packed_bed_reactor()

    with pytest.raises(NotImplementedError):
        reactor._non_isobaric_pressure_balance_gas_phase_reaction()

    with pytest.raises(NotImplementedError):
        reactor._homogeneous_solver()

    with pytest.raises(NotImplementedError):
        reactor._heterogeneous_solver()


def test_asignation_kinetics_arguments():
    "Test for not implemented error of reactor's intherface methods."

    def kinetic1(conc, temperature):
        pass

    def kinetic2(conc, temperature):
        pass

    assert kinetic1([0.25, 0.25, 0.25, 0.25], 298.15) is None
    assert kinetic2([0.25, 0.25, 0.25, 0.25], 298.15) is None

    class SpecificReactor(rd.ReactorBase):
        def __init__(
            self,
            mix,
            list_of_reactions,
            stoichiometry,
            kinetic_argument,
        ):
            self._kinetics = rd.Kinetics(
                list_of_reactions=list_of_reactions,
                mix=mix,
                stoichiometry=stoichiometry,
                kinetic_argument=kinetic_argument,
            )

        @classmethod
        def from_isothermic_isobaric(cls) -> None:
            pass

        @classmethod
        def from_isothermic_noisobaric(cls) -> None:
            pass

        @classmethod
        def from_adiabatic_isobaric(cls) -> None:
            pass

        @classmethod
        def from_adiabatic_noisobaric(cls) -> None:
            pass

        @classmethod
        def from_noisothermic_isobaric(cls) -> None:
            pass

        @classmethod
        def from_noisothermic_noisobaric(cls) -> None:
            pass

        def _set_catalyst_operation(self):
            pass

        def _set_thermal_operation(self):
            pass

        def _set_pressure_operation(self):
            pass

        def _grid_builder(self) -> None:
            pass

        def _border_cond_and_initial_guesses(self) -> None:
            pass

        def _mass_balance(self) -> None:
            pass

        def _energy_balance(self) -> None:
            pass

        def _pressure_balance(self) -> None:
            pass

        def _refrigerant_energy_balance(self) -> None:
            pass

        def simulate(self) -> None:
            pass

        def _homogeneous_mass_balance(self) -> None:
            pass

        def _heterogeneous_mass_balance(self) -> None:
            pass

        def _isothermic_energy_balance(self) -> None:
            pass

        def _homogeneous_adiabatic_energy_balance(self) -> None:
            pass

        def _heterogeneous_adiabatic_energy_balance(self) -> None:
            pass

        def _homogeneous_non_isothermic_energy_balance(self) -> None:
            pass

        def _heterogeneous_non_isothermic_energy_balance(self) -> None:
            pass

        def _isobaric_pressure_balance(self) -> None:
            pass

        def _non_isobaric_pressure_balance_packed_bed_reactor(self) -> None:
            pass

        def _non_isobaric_pressure_balance_gas_phase_reaction(self) -> None:
            pass

        def _homogeneous_solver(self) -> None:
            pass

        def _heterogeneous_solver(self) -> None:
            pass

    assert SpecificReactor.from_isothermic_isobaric() is None
    assert SpecificReactor.from_isothermic_noisobaric() is None
    assert SpecificReactor.from_adiabatic_isobaric() is None
    assert SpecificReactor.from_adiabatic_noisobaric() is None
    assert SpecificReactor.from_noisothermic_isobaric() is None
    assert SpecificReactor.from_noisothermic_noisobaric() is None

    a = rd.Substance(name="SubstanceA")
    b = rd.Substance(name="SubstanceB")

    c = rd.Substance(name="SubstanceC")
    d = rd.Substance(name="SubstanceD")

    mixture1 = rd.mix.IdealSolution(A=a, B=b)
    mixture2 = rd.mix.IdealSolution(C=c, D=d)

    reactor = SpecificReactor(
        mix=mixture1,
        list_of_reactions=[kinetic1],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
    )

    reaction2 = rd.Kinetics(
        mix=mixture2,
        list_of_reactions=[kinetic2],
        stoichiometry=[-2, 2],
        kinetic_argument="partial_pressure",
    )

    # Test for normal calling:

    assert reactor.mix == reactor.kinetics.mix

    assert reactor.list_of_reactions == [kinetic1]
    assert reactor.list_of_reactions == reactor.kinetics.list_of_reactions

    assert np.allclose(reactor.stoichiometry, [-1, 1], atol=1e-30)
    assert np.allclose(
        reactor.stoichiometry, reactor.kinetics.stoichiometry, atol=1e-30
    )

    assert reactor.kinetic_argument == "concentration"
    assert reactor.kinetic_argument == reactor.kinetics.kinetic_argument

    # Changing the kinetics

    reactor.kinetics = reaction2

    assert reactor.mix == reactor.kinetics.mix
    assert reactor.mix == mixture2

    assert reactor.list_of_reactions == [kinetic2]
    assert reactor.list_of_reactions == reactor.kinetics.list_of_reactions

    assert np.allclose(reactor.stoichiometry, [-2, 2], atol=1e-30)
    assert np.allclose(
        reactor.stoichiometry, reactor.kinetics.stoichiometry, atol=1e-30
    )

    assert reactor.kinetic_argument == "partial_pressure"
    assert reactor.kinetic_argument == reactor.kinetics.kinetic_argument

    # change a particular argument

    reactor.stoichiometry = [-1, -1]
    np.allclose(reactor.kinetics.stoichiometry, [-1, -1], atol=1e-30)

    reactor.kinetic_argument = "concentration"
    assert reactor.kinetics.kinetic_argument == "concentration"

    reactor.list_of_reactions = [kinetic1]
    assert reactor.kinetics.list_of_reactions == [kinetic1]

    reactor.mix = mixture1
    assert reactor.kinetics.mix == mixture1

    # making it explode
    with pytest.raises(ValueError):
        reactor.kinetics = 2

    with pytest.raises(ValueError):
        reactor.kinetics = "Messi"

    with pytest.raises(ValueError):
        reactor.kinetics = [1, 2, 3, 4]
