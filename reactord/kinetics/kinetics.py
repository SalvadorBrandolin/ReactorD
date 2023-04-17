from IPython.display import display

from typing import Callable, List

import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic


class Kinetics:
    def __init__(
        self,
        mix: AbstractMix,
        reactive_system: dict,
        rate_argument: str = "concentration",
    ) -> None:
        self.mix = mix
        self.r_argument = rate_argument

        self.r_names: List[str] = list(reactive_system.keys())
        self.r_dicts: List[dict] = [
            reactive_system[name] for name in self.r_names
        ]
        self.r_eqs: List[Symbolic] = [rdict["eq"] for rdict in self.r_dicts]
        self.r_rates: List[Callable] = [
            rdict["rates"] for rdict in self.r_dicts
        ]
        self.r_dh: List[float] = [rdict["DH"] for rdict in self.r_dicts]
        
    def _build_stoichiometry_matrix(self):
        all_names = self.mix.names
        
    def __len__(self):
        return len(self.r_names)

    def __repr__(self):
        ...


d = {"esterificacion": {"eq": asd, "rate": asd, "DH": asd}}
