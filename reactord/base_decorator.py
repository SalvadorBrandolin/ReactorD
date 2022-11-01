from reactord.kinetics import Kinetics
from reactord.reactorbase import ReactorBase


class BaseDecorator(ReactorBase):
    """
    Decorators interface
    """

    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor: ReactorBase = reactor

    @property
    def reactor(self) -> ReactorBase:
        return self._reactor

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
