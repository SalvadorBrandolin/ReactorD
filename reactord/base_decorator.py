from reactord.reactorbase import ReactorBase


class BaseDecorator(ReactorBase):
    """
    Decorators interface
    """
    
    _reactor: ReactorBase = None

    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor = reactor

    @property
    def reactor(self) -> ReactorBase:
        return self._reactor

    # ==================================================================
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self) -> None:
        return self._reactor._grid_builder

    def _border_condition_builder(self) -> None:
        return self._reactor._border_condition_builder

    def _initial_guess_builder(self) -> None:
        return self._reactor._initial_guess_builder

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self) -> None:
        return self._reactor._catalyst_mass_balance

    def _catalyst_energy_balance(self) -> None:
        return self._reactor._catalyst_energy_balance

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self) -> None:
        return self._reactor._mass_balance

    def _reactor_energy_balance(self) -> None:
        return self._reactor._reactor_energy_balance

    def _pressure_balance(self) -> None:
        return self._reactor._pressure_balance

    def _refrigerant_energy_balance(self) -> None:
        return self._reactor._refrigerant_energy_balance

    def simulate(self) -> None:
        return self._reactor.simulate