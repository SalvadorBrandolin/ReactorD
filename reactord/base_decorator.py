import numpy as np

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
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, *args, **kargs) -> None:
        return self._reactor._grid_builder(*args, **kargs)

    def _border_condition_builder(self, *args, **kargs) -> None:
        return self._reactor._border_condition_builder(*args, **kargs)

    def _initial_guess_builder(self, *args, **kargs) -> None:
        return self._reactor._initial_guess_builder(*args, **kargs)

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self, *args, **kargs) -> None:
        return self._reactor._catalyst_mass_balance(*args, **kargs)

    def _catalyst_energy_balance(self, *args, **kargs) -> None:
        return self._reactor._catalyst_energy_balance(*args, **kargs)

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self, *args, **kargs) -> None:
        return self._reactor._mass_balance(*args, **kargs)

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        return self._reactor._reactor_energy_balance(*args, **kargs)

    def _pressure_balance(self, *args, **kargs) -> None:
        return self._reactor._pressure_balance(*args, **kargs)

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        return self._reactor._refrigerant_energy_balance(*args, **kargs)

    def simulate(self, *args, **kargs) -> None:
        return self._reactor.simulate(*args, **kargs)

    # ==================================================================
    # Dunders
    # ==================================================================

    def __getattr__(self, name):
        if not(name in self.__dict__):
            return getattr(self._reactor, name)
        else:
            return object.__getattribute__(self, name)

    def __dir__(self):
        # Append dirs
        list_with_duplicates = dir(self._reactor) + object.__dir__(self)

        # Remove duplicates
        array_without_duplicates = np.unique(list_with_duplicates)

        return list(array_without_duplicates)
         