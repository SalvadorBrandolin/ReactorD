from typing import Callable, List, Union

from IPython.display import Math, display

import numpy as np
from numpy.typing import NDArray

from reactord.kinetic.argument import CompositionalArgument
from reactord.mix.abstract_mix import AbstractMix
from reactord.substance.symbolic import Symbolic

from sympy import latex

from .matrix_builder import stoichiometry_matrix_builder
from .reaction_enthalpy import dh_not_specified, dh_specified


class Kinetic:
    def __init__(
        self,
        mix: AbstractMix,
        reactions: dict,
        kinetic_constants: dict,
        rates_argument: str = "concentration",
    ) -> None:
        # In parameters
        self.mix = mix
        self.r_argument = rates_argument
        self.r_names: List[str] = list(reactions.keys())
        self.r_dics: List[dict] = [reactions[name] for name in self.r_names]
        self.r_eqs: List[Symbolic] = [rdic.get("eq") for rdic in self.r_dics]
        self.r_rates: List[Callable] = [rd.get("rate") for rd in self.r_dics]
        self.r_dh: List[float] = np.array(
            [rdict.get("DH") for rdict in self.r_dics]
        )
        self.kinetic_constants = kinetic_constants

        # Argument evaluation functions
        if self.r_argument == "concentration":
            self._arg_func = self.mix.concentrations
        elif self.r_argument == "partial pressure":
            self._arg_func = self.mix.partial_pressures

        # Build stochiometry matrix from the substances' algebraic expression
        self.stoichiometry = stoichiometry_matrix_builder(self.mix, self.r_eqs)
        
        # Check if all dh are specified, else raise error.
        if any(self.r_dh):
            if not all(self.r_dh):
                raise NotImplementedError(
                    "All reactions must have a reaction enthalpy"
                    " specification or no reaction must have a reaction"
                    " enthalpy specification."
                )
            else:
                self.std_reaction_enthalpies = np.array([])
                self._enthalpy_func = dh_specified
        else:
            hf = self.mix.get_formation_enthalpies()
            self.std_reaction_enthalpies = np.matmul(self.stoichiometry, hf)
            self._enthalpy_func = dh_not_specified

        # Kinetic compositional argument object
        self.comp_argument = CompositionalArgument(self.mix.names)

        # Reaction rates state
        self.r_rates = np.array([])

    def evaluate(
        self,
        mole_fractions: Union[NDArray, float],
        temperature: Union[NDArray, float],
        pressure: Union[NDArray, float],
    ) -> np.ndarray:
        self.comp_argument.values = self._arg_func(
            mole_fractions, temperature, pressure
        )
        self.r_rates = np.array(
            [
                rate(self.comp_argument, temperature, self.kinetic_constants)
                for rate in self.r_rates
            ]
        )
        return self.r_rates

    def dhs_evaluate(self, temperature, pressure):
        return self._enthalpy_func(temperature, pressure)

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
