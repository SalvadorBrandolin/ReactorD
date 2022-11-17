from abc import ABCMeta, abstractmethod

from _collections_abc import Callable

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix


class ReactorBase(metaclass=ABCMeta):
    """Abstract class interface for each reactor in ReactorD."""

    _name: str = ""
    _kinetics: Kinetics = None
    _mix: AbstractMix = None
    _list_of_reactions: list[Callable] = []
    _stoichiometry: list = []
    _kinetic_argument: str = ""

    _catalyst_operation: str = ""
    _thermal_operation: str = ""
    _pressure_operation: str = ""

    _mass_balance_func: Callable = None
    _temperature_balance_func: Callable = None
    _pressure_balance_func: Callable = None

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list[Callable],
        stoichiometry: list,
        kinetic_argument: str,
        **options,
    ) -> None:

        self.options = options
        self.options["_not_reaction_enthalpies"] = True

        self._kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            options=self.options,
        )

    # ==================================================================
    # Common parameters for all reactors.
    # ==================================================================

    @property
    def kinetics(self):
        """Kinetics class instantiation

        Returns
        -------
        Kinetics
            Kinetics class instantiation.
        """
        return self._kinetics

    @kinetics.setter
    def kinetics(self, new_kinetics):
        """Method to asign a new instantiation of the reactor kinetics.
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
    def mix(self):
        """Mixture object that is stored only on the Kinetics object to
        prevent multiple mixture objects on the same reactor.

        Returns
        -------
        AbstractMix
            Mixture object.
        """
        return self._kinetics.mix

    @mix.setter
    def mix(self, new_mix: AbstractMix):
        """Method to assign a new mixture object to the reactor. The
        method replaces the mixture inside the kinetic object.

        Parameters
        ----------
        new_mix : AbastractMix
            Mixture object.
        """
        self._kinetics.mix = new_mix

    @property
    def list_of_reactions(self):
        """List that contains the functions to eval each reaction.

        Explain more: # TODO

        Returns
        -------
        list[Callable]
            List that contains the functions to eval each reaction.
        """
        return self._kinetics.list_of_reactions

    @list_of_reactions.setter
    def list_of_reactions(self, new_list_of_reactions: list[Callable]):
        """Method to assign a new list of reactions to the reactor. The
        method replaces the list of reactions inside the kinetic object.

        Explain more: # TODO

        Parameters
        ----------
        new_list_of_reactions : list[Callable]
            List that contains the functions to eval each reaction.
        """

        self._kinetics.list_of_reactions = new_list_of_reactions

    @property
    def stoichiometry(self):
        """Stoichimetry matrix.

        Explain more: # TODO

        Returns
        -------
        numpy.ndarray[float]
            numpy.ndarray that contains Stoichimetry matrix.
        """

        return self._kinetics.stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        """Stoichimetry matrix replaces on the Kinetics object.

        Explain more: # TODO

        Parameters
        ----------
        new_stoichiometry : list or numpy.ndarray[float]
            Stoichimetry matrix.
        """

        self._kinetics.stoichiometry = new_stoichiometry

    @property
    def kinetic_argument(self):
        """Argument to eval the kinetics function.

        Returns
        -------
        str
            Argument to eval the kinetic functions. Alternatives:
            "concentration", "partial_pressure".
        """

        return self._kinetics.kinetic_argument

    @kinetic_argument.setter
    def kinetic_argument(self, new_kinetics_argument: str):
        """Argument to eval the kinetics function replaced in the
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
    # Abastract methods
    # Settings for mass, energy and pressure balance.
    # ==================================================================

    @abstractmethod
    def set_mass_balance_data(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the mass balance in the reactor's bulk
        phase as attributes of the reactor. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def set_isothermic_operation(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the isothermic energy balance in the
        reactor's bulk. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def set_adiabatic_operation(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the adiabatic energy balance in the
        reactor's bulk. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def set_non_isothermic_operation(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the non isothermic energy balance in the
        reactor's bulk. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def set_isobaric_operation(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the isobaric pressure balance in the
        reactor's bulk phase. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def set_non_isobaric_operation(self):
        """Method that recieves and instantiates the neccesary
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
        """

        Method to build the grid of independent variables.
        Recieves lower and upper boundaries for each independent
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
        """Method that evals and returns the evaluated reactor's bulk
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
        """Method that evals and returns the evaluated reactor's bulk
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
        """Method that evals and returns the evaluated reactor's bulk
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
        """Method that evals and returns the evaluated refrigerant
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
        """Simulates the reactor given the mass_balance_data,
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
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_mass_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics energy balances
    # ==================================================================

    @abstractmethod
    def _isothermic_energy_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _homogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _homogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics pressure balances
    # ==================================================================
    @abstractmethod
    def _isobaric_pressure_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _non_isobaric_pressure_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Specifics solvers
    # ==================================================================
    @abstractmethod
    def _homogeneous_solver(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _heterogeneous_solver(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")
