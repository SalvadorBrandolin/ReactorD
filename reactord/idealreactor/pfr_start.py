from collections.abc import Callable

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase
from reactord.idealreactor.pfr_decorators import StationaryPFR


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

        # Reactor adjectives

        self._adjetives_dictionary = {
            "reactor_type": "Piston flow reactor (PFR)",
            "time_operation": "",
            "catalytic_operation": "",
            "thermal_operation": "",
            "pressure_operation": "",
        }

        # Kinetic set
        self.kinetics: Kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetics_argument=kinetics_argument,
            options=options,
        )

        # Reactors Atributes instantiation
        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

    # ==================================================================
    # Configuration methods
    # ==================================================================

    def set_stationary(
        self, inlet_molar_fluxes: list[float], outlet_molar_fluxes: list[float]
    ):

        self._adjetives["time_operation"] = "stationary"
        
        return StationaryPFR(self)

    def set_non_stationary(self):
        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self) -> None:
        raise NotImplementedError("Text missing")

    def _border_condition_builder(self) -> None:
        raise NotImplementedError("Text missing")

    def _initial_guess_builder(self) -> None:
        raise NotImplementedError("Text missing")

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self) -> None:
        raise NotImplementedError("Text missing")

    def _catalyst_energy_balance(self) -> None:
        raise NotImplementedError("Text missing")

    def _refrigerant_energy_balance(self) -> None:
        raise NotImplementedError("Text missing")

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self) -> None:
        raise NotImplementedError("Text missing")

    def _reactor_energy_balance(self) -> None:
        raise NotImplementedError("Text missing")

    def _pressure_balance(self) -> None:
        raise NotImplementedError("Text missing")

    def simulate(self) -> None:
        raise NotImplementedError("Text missing")
