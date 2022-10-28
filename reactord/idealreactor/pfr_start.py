from collections.abc import Callable
import numpy as np
from scipy.integrate import solve_bvp

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase
from reactord.utils import vectorize


class PFR(ReactorBase):
    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list[Callable],
        stoichiometry: list,
        reactor_dims_minmax: list[float],
        transversal_area: float,
        kinetic_argument: str = "concentration",
    ) -> None:

        # 
        self.kinetic: Kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
        )

        # Reactors Atributes instantiation
        self.reactor_dims_minmax = reactor_dims_minmax
        self.transversal_area = transversal_area

    # ==================================================================
    # Configuration methods
    # ==================================================================
    
    
    
    # ==================================================================
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self) -> None:
        pass

    def _border_condition_builder(self) -> None:
        pass

    def _initial_guess_builder(self) -> None:
        pass

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self) -> None:
        pass

    def _catalyst_energy_balance(self) -> None:
        pass

    def _refrigerant_energy_balance(self) -> None:
        pass

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self) -> None:
        pass

    def _reactor_energy_balance(self) -> None:
        pass

    def _pressure_balance(self) -> None:
        pass

    def simulate(self) -> None:
        pass
