from reactord.kinetics import Kinetics
from reactord.reactorbase import ReactorBase


class BaseDecorator(ReactorBase):
    """
    Decorators interface
    """

    _reactor: ReactorBase = None

    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor: ReactorBase = reactor

    @property
    def reactor(self) -> ReactorBase:
        return self._reactor

    # ==================================================================
    # Specifics for decorators
    # Created to generate shortcus to the decorators attributes
    # ==================================================================

    @property
    def dimensions_settings(self) -> ReactorBase:
        raise NotImplementedError(
            "Time settings are not defined yet. Use '.set_stationary()', or"
            " '.set_non_stationary()'. Check documentation about the reactors"
            " configuration."
        )

    @property
    def time_settings(self) -> ReactorBase:
        raise NotImplementedError(
            "Time settings are not defined yet. Use '.set_stationary()', or"
            " '.set_non_stationary()'. Check documentation about the reactors"
            " configuration."
        )

    @property
    def catalysis_settings(self) -> ReactorBase:
        raise NotImplementedError(
            "Catalysis settings are not defined yet. Use '.set_homogeneous()',"
            " or 'set_heterogeneous()'. Check documentation about the reactors"
            " configuration."
        )

    @property
    def pressure_settings(self) -> ReactorBase:
        raise NotImplementedError(
            "Pressure settings are not defined yet. Use '.set_isobaric()',"
            " or '.set_non_isobaric()'. Check documentation about the reactors"
            " configuration."
        )

    @property
    def thermal_settings(self) -> ReactorBase:
        raise NotImplementedError(
            "Thermal settings are not defined yet. Use '.set_isothermal()',"
            " '.set_adiabatic()' or '.set_non_isothermal()'. Check"
            " documentation about the reactors configuration."
        )

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, *args, **kargs) -> None:
        return self.reactor._grid_builder(*args, **kargs)

    def _border_condition_builder(self, *args, **kargs) -> None:
        return self.reactor._border_condition_builder(*args, **kargs)

    def _initial_guess_builder(self, *args, **kargs) -> None:
        return self.reactor._initial_guess_builder(*args, **kargs)

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self, *args, **kargs) -> None:
        return self.reactor._catalyst_mass_balance(*args, **kargs)

    def _catalyst_energy_balance(self, *args, **kargs) -> None:
        return self.reactor._catalyst_energy_balance(*args, **kargs)

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    @property
    def _settings(self) -> dict:
        return self.reactor._settings

    @_settings.setter
    def _settings(self, updated_settings: dict) -> None:
        self.reactor._settings = dict(updated_settings)

    @property
    def kinetics(self) -> Kinetics:
        return self.reactor.kinetics

    @kinetics.setter
    def kinetics(self, new_kinetics) -> Kinetics:
        self.reactor.kinetics = new_kinetics

    def _mass_balance(self, *args, **kargs) -> None:
        return self.reactor._mass_balance(*args, **kargs)

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        return self.reactor._reactor_energy_balance(*args, **kargs)

    def _pressure_balance(self, *args, **kargs) -> None:
        return self.reactor._pressure_balance(*args, **kargs)

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        return self.reactor._refrigerant_energy_balance(*args, **kargs)

    def simulate(self, *args, **kargs) -> None:
        return self.reactor.simulate(*args, **kargs)

    # ==================================================================
    # Dunders
    # ==================================================================

    def __getattr__(self, name: str):
        return getattr(self.reactor, name)

    def __setattr__(self, name: str, new_value):
        if not (name in self.__dict__):
            return self.reactor.__setattr__(name, new_value)
        else:
            return object.__setattr__(self, name)
