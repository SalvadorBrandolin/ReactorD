"""Reactorbase module."""
from abc import ABCMeta, abstractmethod
from typing import Callable, List

import numpy

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix


class ReactorBase(metaclass=ABCMeta):
    """Abstract class interface for each reactor in ReactorD.

    Parameters
    ----------
    mix : AbstractMix
        Mix object defined with all the substances present in the system
        This object represents properties of mixture of substances
        in the reactor
    list_of_reactions : list[Callable]
        Array that constains kinetic laws for each reaction defined by user
        Laws are add in form of functions like:
        function(composition, temperature)
        where composition is a (number_of_components) dimension array
        that contains the partial pressures [Pa] or the concentrations
        of the substances.
    stoichiometry : list
        array or list containing the stoichiometric coefficients of
        all the substances involved in the reactive system. the
        substances that do not participate in a reaction but are present
        in the reactive system, must have a zero as stoichiometric
        coefficient.

    kinetic_argument : str
        string that indicates on wich concentration unit meassure the
        kinetic rate function are evaluated. Avaliable options:
        'concentration', 'partial_pressure'
    """

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

    # ==================================================================
    # Common parameters for all reactors.
    # ==================================================================

    @property
    def kinetics(self) -> Kinetics:
        """Get Kinetics.

        Return reactor's _kinetic attribute

        Returns
        -------
        Kinetics
            Kinetics class instance.
        """
        return self._kinetics

    @kinetics.setter
    def kinetics(self, new_kinetics: Kinetics) -> None:
        """Set new kinetics object.

        Method to asign a new instantiation of the reactor kinetics.
        Validates that the asigned object is acctualy a Kinetics
        instantiation.

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
                "The assigned value to the reactor's kinetics "
                " attribute must be a Kinetics instance object."
            )

    @property
    def mix(self) -> AbstractMix:
        """Get mix.

        Method to instance Mixture object
        Mixture object is stored only on the Kinetics object to
        prevent multiple mixture objects on the same reactor.

        Returns
        -------
        AbstractMix
            Mixture object.
        """
        return self._kinetics.mix

    @mix.setter
    def mix(self, new_mix: AbstractMix) -> None:
        """Set New Mix object.

        Method to assign a new mixture object to the reactor. The
        method replaces the mixture inside the kinetic object.

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
        Method to store list_of_reactions in Kinetics object

        Returns
        -------
        List[Callable]
            List that contains the functions to eval each reaction.
        """
        return self._kinetics.list_of_reactions

    @list_of_reactions.setter
    def list_of_reactions(self, new_list_of_reactions: List[Callable]):
        """Set new list of reactions.

        Method to assign a new list of reactions to the reactor. The
        method replaces the list of reactions inside the kinetic object.

        Parameters
        ----------
        new_list_of_reactions : List[Callable]
            List that contains the functions to eval each reaction.
        """
        self._kinetics.list_of_reactions = new_list_of_reactions

    @property
    def stoichiometry(self) -> numpy.ndarray:
        """Get stoichiometry.

        Array with stoichimetry information.
        Method to store list_of_reactions in Kinetics object

        Returns
        -------
        numpy.ndarray[float]
            numpy.ndarray that contains Stoichimetry array.
        """
        return self._kinetics.stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry: List[float]) -> None:
        """Set new stoichiometry.

        Stoichimetry array replaces on the Kinetics object.

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
            Argument to eval the kinetic functions. Alternatives:
            "concentration", "partial_pressure".
        """
        return self._kinetics.kinetic_argument

    @kinetic_argument.setter
    def kinetic_argument(self, new_kinetics_argument: str) -> None:
        """Eval kinetic with new argument.

        Argument to eval the kinetics function replaced in the
        Kinetics object.

        Explain better: # TODO

        Parameters
        ----------
        new_kinetics_argument : str
            Argument to eval the kinetic functions. Alternatives:
            "concentration", "partial_pressure".
        """
        self._kinetics.kinetic_argument = new_kinetics_argument

    # ==================================================================
    # Init constructors
    # ==================================================================

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

    # ==================================================================
    # Abastract methods
    # Settings for mass, energy and pressure balance.
    # ==================================================================

    @abstractmethod
    def _set_catalyst_operation(self) -> None:
        """Configure the mass balance settings.

        Method that recieves and instantiates the neccesary
        parameters to solve the mass balance in the reactor's bulk
        phase.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _set_thermal_operation(self) -> None:
        """Configure the energy balance settings.

        Method that receives and instantiates the neccesary
        parameters to solve the energy balance in the
        reactor's bulk phase.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _set_pressure_operation(self):
        """Configure the pressure balance settings.

        Method that receives and instantiates the neccesary
        parameters to solve the non isobaric pressure balance in the
        reactor's bulk phase. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Solvers aditional data needed - general used methods
    # ==================================================================

    @abstractmethod
    def _grid_builder(self) -> None:
        """Build a grid.

        Method to build the grid of independent variables.
        Receives lower and upper boundaries for each independent
        variable and also the number of discretization intervals
        for each defined range.

        Example:

        A discontinuous tank reactor solved from time 0 s to 3600 s,
        using 10 second as time step, should recieve something like:

        time_span = [0, 3600]
        time_step = 10

        The return of the method should seems like:

        return numpy.linspace(time_span[0], time_span[1], time_step)

        More explanation needed: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _border_cond_and_initial_guesses(self) -> None:
        """Construct the solver needs.

        Constructs the border conditions for the differential
        equations that represents the bulk phase behaviour. Also, may
        return intial guesses for the specifics solvers. The format
        of the method's output correponds with the reactor's algebra.

        Explain the output: # TODO

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
        """Eval mass balance.

        Method that evals and returns the evaluated reactor's bulk
        mass balances. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _energy_balance(self) -> None:
        """Eval energy balance.

        Method that evals and returns the evaluated reactor's bulk
        energy balance. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _pressure_balance(self) -> None:
        """Eval pressure balance.

        Method that evals and returns the evaluated reactor's bulk
        pressure balance. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _refrigerant_energy_balance(self) -> None:
        """Eval regrigerant energy balance.

        Method that evals and returns the evaluated refrigerant
        energy balance. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def simulate(self) -> None:
        """Simulate reactor.

        Simulate the reactor given the mass_balance_data,
        energy_balance_data, pressure_balance_data.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics mass balances
    # ==================================================================

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

    # ==================================================================
    # Specifics energy balances
    # ==================================================================

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
        """Eval homogeneous non isothermic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_non_isothermic_energy_balance(self) -> None:
        """Eval heterogeneous non isothermic energy balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics pressure balances
    # ==================================================================
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
        """Eval non isobaric pressure balance.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics solvers
    # ==================================================================
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
