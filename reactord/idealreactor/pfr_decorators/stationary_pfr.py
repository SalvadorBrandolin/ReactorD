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

        self._settings["catalytic_operation"] = "homogenerous"

        return HomogeneousPFR(self)

    def set_heterogeneous(self):
        # TODO
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

    def _grid_builder(self, grid_size: int) -> np.ndarray:
        """Builds the grid for the independent variable "z" (reactor's
        length). The grid goes from the self.reactor_dims_minmax[0] to
        self.reactor_dims_minmax[1] and has a number of elements equal
        to "grid_size".

        Parameters
        ----------
        grid_size : int
            Number of elements of the grid.

        Returns
        -------
        ndarray
            (grid_size, ) dimension ndarray, containing the reactor's
            length grid.
        """
        dim_array = np.linspace(
            self.reactor_dims_minmax[0], self.reactor_dims_minmax[1], grid_size
        )
        return dim_array

    def _border_condition_builder(
        self, ya: np.ndarray, yb: np.ndarray
    ) -> np.ndarray:

        # Checks where the border conditions are given by searching
        # where the values are not np.nan:

        where_not_nan_in = np.invert(np.isnan(self.in_molar_fluxes))
        where_not_nan_out = np.invert(np.isnan(self.in_molar_fluxes))

        inlet = np.argwhere(where_not_nan_in).ravel()
        outlet = np.argwhere(where_not_nan_out).ravel()

        bc_inlet = ya[inlet] - self.in_molar_fluxes[inlet]
        bc_outlet = yb[outlet] - self.out_molar_fluxes[outlet]

        return np.append(bc_inlet, bc_outlet)
