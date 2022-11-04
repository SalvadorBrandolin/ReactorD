import numpy as np

from reactord.decoratorbase import DecoratorBase
from reactord.idealreactor.pfr_decorators import HomogeneousPFR
from reactord.reactorbase import ReactorBase


class StationaryPFR(DecoratorBase):
    def __init__(
        self,
        reactor: ReactorBase,
        inlet_molar_fluxes: list[float],
        outlet_molar_fluxes: list[float],
    ) -> None:

        self._reactor = reactor

        # Specifics attributes
        self.in_molar_fluxes = inlet_molar_fluxes
        self.out_molar_fluxes = outlet_molar_fluxes

        # Validation of fluxes variables

        if np.size(self.in_molar_fluxes) != np.size(self.out_molar_fluxes):
            raise IndexError(
                "set_stationary error: inlet_molar_fluxes list "
                "length must be equal to the outlet_molar_fluxes "
                "list length."
            )

        if np.size(self.in_molar_fluxes) != len(self.kinetic.mix):
            raise IndexError(
                "set_stationary error: inlet_molar_fluxes list "
                "and outlet_molar_fluxes_list length must be "
                "equal to the number of substances in mix "
            )

    # ==================================================================
    # Configuration methods: returns set_catalysis(...(SpecificReactor))
    # ==================================================================

    def set_homogeneous(self) -> HomogeneousPFR:

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "",
                "pressure_operation": "",
            }
        )

        return HomogeneousPFR(self)

    def set_heterogeneous(self):
        """Not Implemented... yet

        Raises
        ------
        NotImplementedError
            method not implemented yet.
        """
        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, grid_size: int):
        """Builds the grid for the independent variable "z" (reactor's length).
        The grid goes from the self.reactor_dims_minmax[0] to
        self.reactor_dims_minmax[1].

        Parameters
        ----------
        grid_size : int
            Number of elements of the grid.

        Returns
        -------
        ndarray
            (grid_size, ) dimension ndarray, containing the reactor's length
            grid.
        """
        dim_array = np.linspace(
            self.reactor_dims_minmax[0], self.reactor_dims_minmax[1], grid_size
        )
        return dim_array

    def _border_condition_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._border_condition_builder(
            *args, **kargs
        )

    def _initial_guess_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._initial_guess_builder(*args, **kargs)

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_mass_balance(*args, **kargs)

    def _catalyst_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_energy_balance(*args, **kargs)

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._mass_balance(*args, **kargs)

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._reactor_energy_balance(*args, **kargs)

    def _pressure_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._pressure_balance(*args, **kargs)

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._refrigerant_energy_balance(
            *args, **kargs
        )

    def simulate(self, *args, **kargs) -> None:
        return self._decorated_reactor.simulate(*args, **kargs)