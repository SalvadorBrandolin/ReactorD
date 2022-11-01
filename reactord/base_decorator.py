from reactord.kinetics import Kinetics
from reactord.reactorbase import ReactorBase


class BaseDecorator(ReactorBase):
    """
    Decorators interface
    """

    _reactor: ReactorBase = None
    _kinetics: Kinetics = None

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
    def kinetics(self) -> Kinetics:
        return self._kinetics

    @property
    def time_settings(self):
        raise NotImplementedError()

    @property
    def catalysis_settings(self):
        raise NotImplementedError()

    @property
    def pressure_settings(self):
        raise NotImplementedError()

    @property
    def thermal_settings(self):
        raise NotImplementedError()

    # ==================================================================
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

    def __getattr__(self, name):
        if not(name in self.__dict__):
            return getattr(self._reactor, name)
        else:
            return object.__getattribute__(self, name)
         