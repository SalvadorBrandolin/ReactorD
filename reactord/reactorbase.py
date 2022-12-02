"""Reactorbase module.

Common interface for reactors.
"""
from abc import ABCMeta, abstractmethod
from typing import Callable, List

import numpy

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix


class ReactorBase(metaclass=ABCMeta):
    """Abstract class interface for each reactor in ReactorD."""

    _name: str = ""
    _kinetics: Kinetics = None
    _mix: AbstractMix = None
    _list_of_reactions: List[Callable] = []
    _stoichiometry: list = []
    _kinetic_argument: str = ""

    _catalyst_operation: str = ""
    _thermal_operation: str = ""
    _pressure_operation: str = ""

    _mass_balance_func: Callable = None
    _energy_balance_func: Callable = None
    _pressure_balance_func: Callable = None
    _solver_func: Callable = None

    # =========================================================================
    # Common parameters for all reactors.
    # =========================================================================

    @property
    def kinetics(self) -> Kinetics:
        """Return reactor's Kinetic object.

        Return reactor's _kinetic attribute.

        Returns
        -------
        Kinetics
            Kinetics class instance.
        """
        return self._kinetics

    @kinetics.setter
    def kinetics(self, new_kinetics: Kinetics) -> None:
        """Set new kinetics object.

        Method to assign a new instance of the reactor kinetics. Validates that
        the assigned object is a Kinetics instance.

        Parameters
        ----------
        new_kinetics : Kinetics
            Kinetics class instantiation.

        Raises
        ------
        ValueError
            The assigned value is not a Kinetics instantiation object.
        """
        if isinstance(new_kinetics, Kinetics):
            self._kinetics = new_kinetics
        else:
            raise ValueError(
                "The assigned value to the reactor's kinetics"
                " attribute must be a Kinetics instance object."
            )

    @property
    def mix(self) -> AbstractMix:
        """Get mix.

        Method to instance Mixture object. The mixture object is stored only
        on the Kinetics object to prevent multiple mixture objects on the same
        reactor.

        Returns
        -------
        AbstractMix
            Mixture object.
        """
        return self._kinetics.mix

    @mix.setter
    def mix(self, new_mix: AbstractMix) -> None:
        """Set New Mix object.

        Method to assign a new mixture object to the reactor. The method
        replaces the mixture inside the kinetic object.

        Parameters
        ----------
        new_mix : AbastractMix
            Mixture object.

        Returns
        -------
        AbstractMix
            Mixture object.
        """
        self._kinetics.mix = new_mix

    @property
    def list_of_reactions(self) -> List[Callable]:
        """Get list of reactions.

        List that contains the functions to eval each law kinetic reaction.
        Method to store list_of_reactions in Kinetics object.

        Returns
        -------
        List[Callable]
            List that contains the functions to eval each reaction.
        """
        return self._kinetics.list_of_reactions

    @list_of_reactions.setter
    def list_of_reactions(self, new_list_of_reactions: List[Callable]):
        """Set a new list of reactions.

        Method to assign a new list of reactions to the reactor. The method
        replaces the list of reactions inside the kinetic object.

        Parameters
        ----------
        new_list_of_reactions : List[Callable]
            List that contains the functions to eval each reaction.
        """
        self._kinetics.list_of_reactions = new_list_of_reactions

    @property
    def stoichiometry(self) -> numpy.ndarray:
        """Get stoichiometry.

        Numpy.ndrray with stoichiometry information. Method to store
        list_of_reactions in Kinetics object.

        Returns
        -------
        numpy.ndarray[float]
            numpy.ndarray that contains the Stoichiometry array.
        """
        return self._kinetics.stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry: List[float]) -> None:
        """Set new stoichiometry.

        The new_stoichimetry array replaces the reactor's kinetic
        stoichiometry.

        Parameters
        ----------
        new_stoichiometry : list or numpy.ndarray[float]
            Stoichimetry matrix.
        """
        self._kinetics.stoichiometry = new_stoichiometry

    @property
    def kinetic_argument(self) -> str:
        """Argument to eval the kinetics function.

        Returns
        -------
        str
            Argument to eval the kinetic laws. Alternatives:
            "concentration", "partial_pressure".
        """
        return self._kinetics.kinetic_argument

    @kinetic_argument.setter
    def kinetic_argument(self, new_kinetics_argument: str) -> None:
        """Replace kinetic argument with a new argument.

        Argument to eval the kinetics function replaced in the Kinetics object.

        Parameters
        ----------
        new_kinetics_argument : str
            Argument to eval the kinetic functions. Alternatives:
            "concentration", "partial_pressure".
        """
        self._kinetics.kinetic_argument = new_kinetics_argument

    # =========================================================================
    # Init constructors
    # =========================================================================

    @classmethod
    @abstractmethod
    def set_isothermic_isobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @classmethod
    @abstractmethod
    def set_isothermic_noisobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @classmethod
    @abstractmethod
    def set_adiabatic_isobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @classmethod
    @abstractmethod
    def set_adiabatic_noisobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @classmethod
    @abstractmethod
    def set_noisothermic_isobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    @classmethod
    @abstractmethod
    def set_noisothermic_noisobaric(cls) -> None:
        """Abstract method not implemented.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented")

    # =========================================================================
    # Abastract methods
    # Settings for mass, energy and pressure balance.
    # =========================================================================

    @abstractmethod
    def _set_catalyst_operation(self) -> None:
        """Configure the mass balance settings.

        Method that receives and instantiates the necessary parameters to solve
        the mass balance in the reactor's bulk phase.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _set_thermal_operation(self) -> None:
        """Configure the energy balance settings.

        Method that receives and instantiates the necessary parameters to solve
        the energy balance in the reactor's bulk phase.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _set_pressure_operation(self) -> None:
        """Configure the pressure balance settings.

        Method that receives and instantiates the necessary parameters to solve
        the non-isobaric pressure balance in the reactor's bulk phase.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # =========================================================================
    # Solvers additional data needed - generally used methods
    # =========================================================================

    @abstractmethod
    def _grid_builder(self) -> None:
        """Build the independent variables a grid.

        Method to build the grid of independent variables. Receives lower and
        upper boundaries for each independent variable and also the number of
        discretization intervals for each defined range.

        Example:

        A discontinuous tank reactor solved from time 0 s to 3600 s, using 10
        seconds as a time step should receive something like:

        time_span = [0, 3600]
        time_step = 10

        The return of the method should seem like this:

        return numpy.linspace(time_span[0], time_span[1], time_step)

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _border_cond_and_initial_guesses(self) -> None:
        """Construct the solver's border conditions and initial guesses.

        Constructs the border conditions for the differential equations that
        represent the bulk phase behavior. Also, may return initial guesses for
        the specific solvers. The format of the method's output corresponds
        with the reactor's algebra.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Balances
    # ==================================================================

    @abstractmethod
    def _mass_balance(self) -> None:
        """Eval the reactor's mass balance.

        Method that evals and returns the evaluated reactor's bulk mass
        balances. The format of the method's returns corresponds to the
        specific solver's needs.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _energy_balance(self) -> None:
        """Eval reactor's energy balance.

        Method that evals and returns the evaluated reactor's bulk energy
        balance. The format of the method's returns corresponds to the specific
        solver's needs.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _pressure_balance(self) -> None:
        """Eval reactor's pressure balance.

        Method that evals and returns the evaluated reactor's bulk pressure
        balance. The format of the method's returns corresponds to the specific
        solver's needs.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _refrigerant_energy_balance(self) -> None:
        """Eval regrigerant energy balance.

        Method that evals and returns the evaluated refrigerant energy balance.
        The format of the method's returns corresponds to the specific solver's
        needs.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def simulate(self) -> None:
        """Simulate reactor.

        Simulate the reactor given the mass_balance_data, energy_balance_data,
        pressure_balance_data.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # =========================================================================
    # Specifics mass balances
    # =========================================================================

    @abstractmethod
    def _homogeneous_mass_balance(self) -> None:
        """Eval homogeneous mass balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_mass_balance(self) -> None:
        """Eval heterogeneous mass balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # =========================================================================
    # Specifics energy balances
    # =========================================================================

    @abstractmethod
    def _isothermic_energy_balance(self) -> None:
        """Eval isothermic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _homogeneous_adiabatic_energy_balance(self) -> None:
        """Eval homogeneous adiabatic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_adiabatic_energy_balance(self) -> None:
        """Eval heterogeneous adiabatic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _homogeneous_non_isothermic_energy_balance(self) -> None:
        """Eval homogeneous non-isothermic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_non_isothermic_energy_balance(self) -> None:
        """Eval heterogeneous non-isothermic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # =========================================================================
    # Specifics pressure balances
    # =========================================================================
    @abstractmethod
    def _isobaric_pressure_balance(self) -> None:
        """Eval isobaric pressure balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _non_isobaric_pressure_balance(self) -> None:
        """Eval non-isobaric pressure balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # =========================================================================
    # Specifics solvers
    # =========================================================================
    @abstractmethod
    def _homogeneous_solver(self) -> None:
        """Solve homogeneous reactor.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_solver(self) -> None:
        """Solve heterogeneous reactor.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")
