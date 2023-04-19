from IPython.display import display

from typing import Callable, List

import numpy as np

from .matrix_builder import stoichiometry_matrix_builder

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic


class Kinetics:
    def __init__(
        self,
        mix: AbstractMix,
        reactions: dict,
        rates_argument: str = "concentration",
    ) -> None:
        self.mix = mix
        self.r_argument = rates_argument

        self.r_names: List[str] = list(reactions.keys())
        self.r_dics: List[dict] = [reactions[name] for name in self.r_names]
        self.r_eqs: List[Symbolic] = [rdic.get("eq") for rdic in self.r_dics]
        self.r_rates: List[Callable] = [rd.get("rate") for rd in self.r_dics]
        self.r_dh: List[float] = [rdict.get("DH") for rdict in self.r_dics]

        self.stoichiometry = stoichiometry_matrix_builder(self.mix, self.r_eqs)

    def _build_stoichiometry_matrix(self):
        all_names = self.mix.names

    def __len__(self):
        return len(self.r_names)

    def __repr__(self):
        ...
