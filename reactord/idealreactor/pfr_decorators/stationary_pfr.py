import numpy as np

from reactord.base_decorator import BaseDecorator
from reactord.idealreactor import PFR
from reactord.kinetics import Kinetics


class StationaryPFR(BaseDecorator):
    def __init__(
        self,
        reactor: PFR,
        inlet_molar_fluxes: list[float],
        outlet_molar_fluxes: list[float],
    ) -> None:

        self._reactor = reactor
        self.inlet_molar_fluxes = inlet_molar_fluxes
        self.inlet_molar_fluxes = outlet_molar_fluxes

    # ==================================================================
    # Configuration methods: returns set_catalysis(set_time(PFR))
    # ==================================================================
    
    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, grid_size: int):
        dim_array = np.linspace(
            self.reactor_dims_minmax[0], self.reactor_dims_minmax[1], grid_size
        )
        return dim_array

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
