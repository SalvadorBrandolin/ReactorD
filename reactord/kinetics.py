"""Kinetics module."""
from typing import List

import numpy as np

from .mix.abstract_mix import AbstractMix


class Kinetics:
    """Kinetic object builder.

    Parameters
    ----------
    mix : AbstractMix
        Mixture object.
    list_of_reactions : List[Callable]
        List of functions that evaluate the reaction rates, each one
        defined by the user with the following format:
        callable(concentration_unit: list[float], temperature: float
        ) -> float.
        Where concentration_unit refers to the units in which the
        arguments of the kinetic laws are expressed, for instance,
        concentrations or partial_pressures.
    stoichiometry : List[float]
        A matrix that represents the stoichiometry of the reactive system.
        Each row represents one reaction contained in the
        list_of_reactions parameter and each column represents a
        substance in the mix parameter. The stoichiometry matrix
        entrances are the stoichiometric coefficients of each substance
        in each reaction.
    kinetic_argument : str
        This argument is used to define how to evaluate the composition
        of the reactive mixture inside the reactor.
        Options:
        'concentration': substance concentration. [mol/m^3]
        'partial_pressure': substance partial pressure. [Pa]

    Attributes
    ----------
    mix : AbstractMix
        Mixture object.
    list_of_reactions : List[Callable]
        List of functions that evaluate the reaction rates, each one
        defined by the user with the following format:
        callable(concentration_unit: list[float], temperature: float
        ) -> float.
        Where concentration_unit refers to the units in which the
        arguments of the kinetic laws are expressed, for instance,
        concentrations or partial_pressures.
    stoichiometry : ndarray[float]
        A matrix that represents the stoichiometry of the reactive system.
        Each row represents one reaction contained in the
        list_of_reactions parameter and each column represents a
        substance in the mix parameter. The stoichiometry matrix
        entrances are the stoichiometric coefficients of each substance
        in each reaction.
    kinetic_argument : str
        This argument is used to define how to evaluate the composition
        of the reactive mixture inside the reactor.
        Options:
        'concentration': substance concentration. [mol/m^3]
        'partial_pressure': substance partial pressure. [Pa]
    user_reaction_enthalpies : ndarray[float]
        User specified reaction enthalpies.
    num_reactions : int
        Number of reactions.
    num_substances : int
        Number of substances in mix.
    """

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list,
        stoichiometry: list,
        kinetic_argument: str = "concentration",
        reaction_enthalpies=None,
    ) -> None:

        self.list_of_reactions = list_of_reactions
        self.mix = mix
        self.kinetic_argument = kinetic_argument.lower()
        self.user_reaction_enthalpies = None
        self._std_reaction_enthalpies = None

        # =====================================================================
        # DATA VALIDATION
        # =====================================================================

        # Is mix an instance of AbstractMix?
        if not (isinstance(self.mix, AbstractMix)):
            raise TypeError(
                "The supplied argument 'mix' must be an instance of mixture "
                "class. See documentation for defining mix objects."
            )

        # Get the number of components and reactions from stoichiometry
        if np.ndim(stoichiometry) == 1:
            self.num_reactions = 1
            self.num_substances = np.shape(stoichiometry)[0]
        else:
            self.num_reactions, self.num_substances = np.shape(stoichiometry)

        # There must be one function to evaluate the reaction rate per reaction
        if self.num_reactions != len(list_of_reactions):
            raise IndexError(
                "'stoichiometry' rows number must be equal to"
                " list_of_reactions' length"
            )

        # Check that mix and stoichiometry have the same number of substances
        if len(mix) != self.num_substances:
            raise IndexError(
                "'stoichiometry' columns number must be equal to substances"
                " number in 'mix' object"
            )

        # Checks whether reaction_enthalpies option in kwargs is correct
        if reaction_enthalpies is not None:
            if len(reaction_enthalpies) != self.num_reactions:
                raise IndexError(
                    "The number of reaction enthalpies in the"
                    " reaction_enthalpies option"
                    f" [{len(reaction_enthalpies)}] must be equal to"
                    " the stoichiometry matrix row number"
                )

            self.user_reaction_enthalpies = np.array(reaction_enthalpies)

        # =====================================================================
        # Set the dimension of the stoichiometry matrix explicitly
        # (needed for single reaction systems)
        # =====================================================================
        self.stoichiometry = np.array(stoichiometry).reshape(
            self.num_reactions, self.num_substances
        )

        # =====================================================================
        # Set kinetics compositional argument
        # =====================================================================
        if self.kinetic_argument == "concentration":
            self._composition_calculator = self.mix.concentrations
        elif self.kinetic_argument == "partial_pressure":
            self._composition_calculator = self.mix.partial_pressures
        else:
            raise ValueError(
                f"{self.kinetic_argument} is not a valid kinetic argument"
            )

    # ==================================================================
    # PUBLIC METHODS
    # ==================================================================
    @property
    def std_reaction_enthalpies(self) -> np.ndarray:
        """Set standard reaction enthalpies.

        Array containing the standard reaction enthalpies of the Kinetic
        object. This attribute must be initializated by the
        std_reaction_enthalpies_init method, by default None.
        """
        return self._std_reaction_enthalpies

    @std_reaction_enthalpies.setter
    def std_reaction_enthalpies(self, new_kinetics: List[float]) -> None:
        raise NotImplementedError(
            "The attribute std_reaction_enthalpies doesn't admit direct"
            "assignation."
        )

    def kinetic_eval(
        self, moles: List[float], temperature: float, pressure: float
    ) -> np.ndarray:
        """Evaluate reaction kinetics.

        Method that evaluates the reaction rates for each reaction and
        each component in the mixture at given concentrations, temperature
        and pressure.

        Parameters
        ----------
        moles : ndarray or list
            moles of each substance
        temperature : float
            Temperature [K]
        pressure : float
            Pressure [Pa]

        Returns
        -------
        rates_i: ndarray of shape (moles,). The rates for each component.
        reaction_rates: ndarray of shape (list_reactions,). The rate for
            each reaction.
        """
        # The partial pressures or concentrations are calculated:
        composition = self._composition_calculator(
            moles, temperature, pressure
        )

        # Rates for each reaction:
        reaction_rates = np.array(
            [
                reaction(composition, temperature)
                for reaction in self.list_of_reactions
            ]
        )

        # Rates for each compound:
        rates_i = np.matmul(reaction_rates, self.stoichiometry)
        return rates_i, reaction_rates

    def std_reaction_enthalpies_init(self) -> None:
        """Calculate the standard reaction enthalpies.

        If reaction enthalpies were not specified in the Kinetic object
        initialization, this method will calculate the standard reaction
        enthalpies from the standard formation enthalpies stored in the pure
        compounds of the mix.

        This method is called only when reaction enthalpies are necessary for
        the reactors' non-isothermic operations. Allowing that no information
        on the standard formation enthalpies and heat capacities of the pure
        substances is necessary when an isothermic operation is performed.
        """
        if self.user_reaction_enthalpies is not None:
            pass
        else:
            self._std_reaction_enthalpies = (
                self._std_reaction_enthalpies_from_formation()
            )

    def reaction_enthalpies(
        self, temperature: float, pressure: float
    ) -> np.ndarray:
        """Evaluate reaction enthalpies of all the reactions involved.

        Parameters
        ----------
        temperature : float
            Temperature [K] at which enthalpies are evaluated.
        pressure : float
            Pressure [Pa] at which enthalpies are evaluated.

        Returns
        -------
        ndarray
            Array containint the reaction enthalpies at the specified
            temperature. [J/mol]
        """
        if self.user_reaction_enthalpies is not None:
            return self.user_reaction_enthalpies

        else:

            formation_correction = self.mix.formation_enthalpies_correction(
                temperature, pressure
            )

            reaction_enthalpies_correction = np.dot(
                self.stoichiometry, formation_correction
            )

            reaction_enthalpies = (
                reaction_enthalpies_correction + self.std_reaction_enthalpies
            )
            return reaction_enthalpies

    # =========================================================================
    # PRIVATE METHODS
    # =========================================================================
    def _std_reaction_enthalpies_from_formation(self) -> None:
        """Calculate standard reaction enthalpies.

        Calculates the standard reaction enthalpies using standard
        formation enthalpies from a mix object.

        Returns
        -------
        ndarray
            Array containing the standard reaction enthalpies. [J/mol/K]
        """
        formation_enthalpies = self.mix._formation_enthalpies_set()
        return np.dot(self.stoichiometry, formation_enthalpies)
