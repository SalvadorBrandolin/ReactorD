"""Kinetics module."""
import numpy as np

from .mix.abstract_mix import AbstractMix
from .utils import vectorize


class Kinetics:
    """Kinetic object builder.

    Parameters
    ----------
    mix : AbstractMix
            Mixture object.
    list_of_reactions : List[Callable]
        List containing functions to eval the reaction rates, each
        defined by the user with the form:
        callable(concentration_unit: list[float], temperature: float
        ) -> float. Where concentration_unit refers to the
        concentration unit of measure that is argument of the
        kinetic law.
    stoichiometry : List[float]
        Stoichiometry matrix of the reactive system. Each row
        represents each reaction contained in the list_of_reactions
        parameter, each column represents each substance in the mix
        parameter. The stoichiometry matrix entrances are the
        stoichiometric coefficients of each substance in each
        reaction.
    kinetic_argument : str
        Kinetic argument used to eval the reaction defined by the
        user. Options:
        'concentration': substance concentration. [mol/m^3]
        'partial_pressure': substance partial pressure. [Pa]
    """

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list,
        stoichiometry: list,
        kinetic_argument: str = "concentration",
        **kwargs,
    ) -> None:

        self.list_of_reactions = list_of_reactions
        self.mix = mix
        self.kinetic_argument = kinetic_argument.lower()
        self.kwargs = kwargs

        # ==============================================================
        # DATA VALIDATION
        # ==============================================================

        # Mix is an instances of AbstractMix?

        if not (isinstance(self.mix, AbstractMix)):
            raise TypeError(
                "The supplied argument 'mix' must be an instance of mixture "
                "class. See documentation for defining mix objects."
            )

        # Get the number of component and reactions from stoichiometry

        if np.ndim(stoichiometry) == 1:
            self.num_reactions = 1
            self.num_substances = np.shape(stoichiometry)[0]
        else:
            self.num_reactions, self.num_substances = np.shape(stoichiometry)

        # Checks that there is a function to eval each reaction

        if self.num_reactions != len(list_of_reactions):
            raise IndexError(
                "'stoichiometry' rows number must be equal to"
                " list_of_reactions' length"
            )

        # Checks mix and stoichiometry has the same number of substances

        if len(mix) != self.num_substances:
            raise IndexError(
                "'stoichiometry' columns number must be equal to substances"
                " number in 'mix' object"
            )

        # Checks if reacion_enthalpies option is correct

        if "reaction_enthalpies" in self.kwargs.keys():

            reaction_enthalpies = self.kwargs.get("reaction_enthalpies")

            if len(reaction_enthalpies) != self.num_reactions:
                raise IndexError(
                    "The number of reaction enthalpies in the"
                    " reaction_enthalpies option"
                    f" [{len(reaction_enthalpies)}] must be equal to"
                    " the stoichiometry matrix row number"
                )

        # ==============================================================
        # Set the dimension of stoichiometry matrix explicitly
        # (needed for single reaction systems)
        # ==============================================================

        self.stoichiometry = np.array(stoichiometry).reshape(
            self.num_reactions, self.num_substances
        )

        # ==============================================================
        # SETS THE KINETICS COMPOSITIONAL ARGUMENTS
        # ==============================================================

        if self.kinetic_argument == "concentration":
            self._composition_calculator = self.mix.concentrations

        elif self.kinetic_argument == "partial_pressure":
            self._composition_calculator = self.mix.partial_pressures

        else:
            raise ValueError(
                f"{self.kinetic_argument} is not a valid kinetic argument"
            )

        # ==============================================================
        # FORMATION AND REACTION ENTHALPIES SET
        # ==============================================================

        self._std_reaction_enthalpies = None

    # ==================================================================
    # PUBLIC METHODS
    # ==================================================================

    @vectorize(signature="(n),(),()->(n),(m)", excluded={0})
    def kinetic_eval(
        self, moles: list, temperature: float, pressure: float
    ) -> np.ndarray:
        """Evaluate kinetic.

        Method that evaluates the reaction rate for each reaction at
        concentration, temperature and pressure given.

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
        ndarray, ndarray

        """
        # The partial pressures or concentrations are calculated:
        composition = self._composition_calculator(
            moles, temperature, pressure
        )

        # Rates for each individual reaction:
        reaction_rates = np.array(
            [
                reaction(composition, temperature)
                for reaction in self.list_of_reactions
            ]
        )
        # Rates for each compound:
        rates_i = np.matmul(reaction_rates, self.stoichiometry)
        return rates_i, reaction_rates

    @property
    def std_reaction_enthalpies(self):
        """Set standar reaction enthalpies."""
        return self._std_reaction_enthalpies

    @std_reaction_enthalpies.setter
    def std_reaction_enthalpies(self):
        raise NotImplementedError(
            "The attribute std_reaction_enthalpies doesn't admits direct"
            "assignation."
        )

    @vectorize(signature="(),()->(m)", excluded={0})
    def reaction_enthalpies(self, temperature, pressure):
        """Eval reacion enthalpies.

        Parameters
        ----------
        temperature : float
            temperature to eval enthalpies
        pressure : float
            pressure to eval enthalpies

        Returns
        -------
        array, attribute
            reaction enthalpies with correction acording to the substance and
            temperature
        """
        if "reaction_enthalpies" in self.kwargs.keys():
            return self.kwargs.get("reaction_enthalpies")

        formation_correction = self.mix.formation_enthalpies_correction(
            temperature, pressure
        )

        reaction_enthalpies_correction = np.dot(
            self.stoichiometry, formation_correction
        )

        return reaction_enthalpies_correction + self.std_reaction_enthalpies

    # ==================================================================
    # PRIVATE METHODS
    # ==================================================================

    def _std_reaction_enthalpies_set(self):
        if "reaction_enthalpies" in self.kwargs.keys():
            pass
        else:
            self._std_reaction_enthalpies = (
                self._std_reaction_enthalpies_from_formation()
            )

    def _std_reaction_enthalpies_from_formation(self):
        """Calculate standar reaction enthalpies.

        Calculates the standard reaction enthalpy from standard
        formation enthalpies defined from mix object.

        Returns
        -------
        ndarray
            Standar reaction enthalpies. [j/mol/K]
        """
        formation_enthalpies = self.mix._formation_enthalpies_set()
        return np.dot(self.stoichiometry, formation_enthalpies)
