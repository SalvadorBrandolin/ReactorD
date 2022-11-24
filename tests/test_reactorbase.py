import numpy as np

import pytest

import reactord as rd


def test_abstract_error_raises():

    with pytest.raises(TypeError):

        def kinetic(conc, temperature):
            return 10

        a = rd.Substance()
        b = rd.Substance()

        mixture = rd.mix.IdealSolution([a, b])

        reactor = rd.ReactorBase(
            mix=mixture,
            list_of_reactions=[kinetic],
            stoichiometry=[-1, 1],
            kinetic_argument="concentration",
        )

        reactor._mass_balance_func()


def test_heritage_methods():

    # ==================================================================
    # Not implemented all abstract methods
    # ==================================================================
    class SpecificReactor(rd.ReactorBase):
        ...

    with pytest.raises(TypeError):
        reactor = SpecificReactor()
        reactor._mass_balance_func()

    # ==================================================================
    # Implemented some abstract methods
    # ==================================================================

    class SpecificReactor(rd.ReactorBase):
        def _grid_builder(self) -> None:
            ...

        def _border_cond_and_initial_guesses(self) -> None:
            ...

    with pytest.raises(TypeError):
        reactor = SpecificReactor()
        reactor._mass_balance_func()


def test_not_implemented_error():
    "Test for not implemented error of reactor's intherphase methods."

    def kinetic(conc, temperature):
        return 10

    class SpecificReactor(rd.ReactorBase):
        def set_mass_balance_data(self):
            return super().set_mass_balance_data()

        def set_isothermic_operation(self):
            return super().set_isothermic_operation()

        def set_adiabatic_operation(self):
            return super().set_adiabatic_operation()

        def set_non_isothermic_operation(self):
            return super().set_non_isothermic_operation()

        def set_isobaric_operation(self):
            return super().set_isobaric_operation()

        def set_non_isobaric_operation(self):
            return super().set_non_isobaric_operation()

        def _grid_builder(self) -> None:
            return super()._grid_builder()

        def _border_cond_and_initial_guesses(self) -> None:
            return super()._border_cond_and_initial_guesses()

        def _mass_balance(self) -> None:
            return super()._mass_balance()

        def _energy_balance(self) -> None:
            return super()._energy_balance()

        def _pressure_balance(self) -> None:
            return super()._pressure_balance()

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

        def _non_isobaric_pressure_balance(self) -> None:
            return super()._non_isobaric_pressure_balance()

        def _homogeneous_solver(self) -> None:
            return super()._homogeneous_solver()

        def _heterogeneous_solver(self) -> None:
            return super()._heterogeneous_solver()

    a = rd.Substance()
    b = rd.Substance()

    mixture = rd.mix.IdealSolution([a, b])

    reactor = SpecificReactor(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
    )

    with pytest.raises(NotImplementedError):
        reactor.set_mass_balance_data()

    with pytest.raises(NotImplementedError):
        reactor.set_isothermic_operation()

    with pytest.raises(NotImplementedError):
        reactor.set_adiabatic_operation()

    with pytest.raises(NotImplementedError):
        reactor.set_non_isothermic_operation()

    with pytest.raises(NotImplementedError):
        reactor.set_isobaric_operation()

    with pytest.raises(NotImplementedError):
        reactor.set_non_isobaric_operation()

    with pytest.raises(NotImplementedError):
        reactor._grid_builder()

    with pytest.raises(NotImplementedError):
        reactor._border_cond_and_initial_guesses()

    with pytest.raises(NotImplementedError):
        reactor._mass_balance()

    with pytest.raises(NotImplementedError):
        reactor._energy_balance()

    with pytest.raises(NotImplementedError):
        reactor._pressure_balance()

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
        reactor._non_isobaric_pressure_balance()

    with pytest.raises(NotImplementedError):
        reactor._homogeneous_solver()

    with pytest.raises(NotImplementedError):
        reactor._heterogeneous_solver()


def test_asignation_kinetics_arguments():
    "Test for not implemented error of reactor's intherphase methods."

    def kinetic1(conc, temperature):
        return 100

    def kinetic2(conc, temperature):
        return 200

    class SpecificReactor(rd.ReactorBase):
        def set_mass_balance_data(self):
            return super().set_mass_balance_data()

        def set_isothermic_operation(self):
            return super().set_isothermic_operation()

        def set_adiabatic_operation(self):
            return super().set_adiabatic_operation()

        def set_non_isothermic_operation(self):
            return super().set_non_isothermic_operation()

        def set_isobaric_operation(self):
            return super().set_isobaric_operation()

        def set_non_isobaric_operation(self):
            return super().set_non_isobaric_operation()

        def _grid_builder(self) -> None:
            return super()._grid_builder()

        def _border_cond_and_initial_guesses(self) -> None:
            return super()._border_cond_and_initial_guesses()

        def _mass_balance(self) -> None:
            return super()._mass_balance()

        def _energy_balance(self) -> None:
            return super()._energy_balance()

        def _pressure_balance(self) -> None:
            return super()._pressure_balance()

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

        def _non_isobaric_pressure_balance(self) -> None:
            return super()._non_isobaric_pressure_balance()

        def _homogeneous_solver(self) -> None:
            return super()._homogeneous_solver()

        def _heterogeneous_solver(self) -> None:
            return super()._heterogeneous_solver()

    a = rd.Substance(name="SubstanceA")
    b = rd.Substance(name="SubstanceB")

    c = rd.Substance(name="SubstanceC")
    d = rd.Substance(name="SubstanceD")

    mixture1 = rd.mix.IdealSolution([a, b])
    mixture2 = rd.mix.IdealSolution([c, d])

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
