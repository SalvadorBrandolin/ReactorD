from IPython.display import display, Math, Markdown

from typing import Callable, List

import numpy as np

from .matrix_builder import stoichiometry_matrix_builder

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic
from reactord.kinetic.argument import KineticArgument

from sympy import latex


class Kinetic:
    def __init__(
        self,
        mix: AbstractMix,
        reactions: dict,
        rates_argument: str = "concentration",
        **kinetic_constants
    ) -> None:
        
        self.mix = mix
        self.r_argument = rates_argument
        self.r_names: List[str] = list(reactions.keys())
        self.r_dics: List[dict] = [reactions[name] for name in self.r_names]
        self.r_eqs: List[Symbolic] = [rdic.get("eq") for rdic in self.r_dics]
        self.r_rates: List[Callable] = [rd.get("rate") for rd in self.r_dics]
        self.r_dh: List[float] = [rdict.get("DH") for rdict in self.r_dics]
        self.kinetic_constants = kinetic_constants

        self.stoichiometry = stoichiometry_matrix_builder(self.mix, self.r_eqs)
        self.arguments = KineticArgument(self.mix.names)

    @property
    def irepr(self):
        for r_name, eq in zip(self.r_names, self.r_eqs):
            ltx = latex(eq._chem_equality).replace("=", r"\rightarrow")
            display(Math(f"{r_name}: {ltx}"))

    def __len__(self):
        return len(self.r_names)
            
    def __repr__(self):
        output = "Mixture's substances: \n"
        
        for name in self.mix.names:
            output = output + f"  * {name} \n"
        
        output = output + "\n"
        output = output + "System's reactions: \n"
        
        for r_name, eq in zip(self.r_names, self.r_eqs):
            ltx = latex(eq._chem_equality).replace("=", r"\rightarrow")
            output = output + f"{r_name}: {ltx} \n"
        
        return output
            
        
