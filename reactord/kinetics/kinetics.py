from typing import Callable, List, Union

import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic


class Kinetics:
    def __init__(self, mix: AbstractMix, reactive_system: List[dict]) -> None:
        self.mix = mix
        self.exprs: List[Symbolic] = [r["eq"] for r in reactive_system]
        self.r_list = List[Callable] = [r["rate"] for r in reactive_system]
        self.dhs = List[Union[None, float]] = [
            r["DH"] for r in reactive_system
        ]

