from collections.abc import Callable

from reactord.idealreactor.pfr_decorators import StationaryPFR
from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase


class PFR(ReactorBase):
    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list[Callable],
        stoichiometry: list,
        kinetics_argument: str,
        reactor_dim_minmax: list[float],
        transversal_area: float,
        **options,
    ) -> None:

        # Kinetic set
        self.kinetics: Kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetics_argument=kinetics_argument,
            options=options,
        )

        # Reactors dimentions
        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

        # Settings
        self._settings = {
            "reactor_type": "Piston flow reactor (PFR)",
            "time_operation": "",
            "catalytic_operation": "",
            "thermal_operation": "",
            "pressure_operation": "",
        }

    # ==================================================================
    # Configuration methods: returns set_time(SpecificReactor)
    # ==================================================================

    def set_stationary(
        self, inlet_molar_fluxes: list[float], outlet_molar_fluxes: list[float]
    ) -> StationaryPFR:

        self._settings["time_operation"] = "stationary"

        return StationaryPFR(self, inlet_molar_fluxes, outlet_molar_fluxes)

    def set_non_stationary(self):
        # TODO
        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _border_condition_builder(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _initial_guess_builder(self,*args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _catalyst_energy_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _catalyst_border_conditions(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _catalyst_initial_guess_builder(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _pressure_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")

    def simulate(self, *args, **kargs) -> None:
        raise NotImplementedError("Text missing")
