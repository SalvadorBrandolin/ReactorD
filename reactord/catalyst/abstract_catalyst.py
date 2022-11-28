from abc import ABCMeta, abstractmethod


class AbstractCatalyst(metaclass=ABCMeta):
    """Abstract class interface for particles in
    heterogeneous reactors in ReactorD."""

    # ==================================================================
    # Abastract methods
    # ODE/PDE particles general used methods
    # ==================================================================

    @abstractmethod
    def _grid_builder_particle(self) -> None:
        """Method to build the grid of independent variables.
        Recieves lower and upper boundaries for each independent
        variable and also the number of discretization intervals
        for each defined range.

        Example:

        A sphere particle solved from 0 to radio,
        using 0.1 as step, should recieve something like:

        dimension_span = [0, 1]
        dimension_step = 0.1

        The return of the method should seems like:

        return numpy.linspace(dimension_span[0], dimension_span[1],
        dimension_step)
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _border_condition_builder_particle(self) -> None:

        """Constructs the border conditions for the differential
        equations that represents the particle phase behaviour.
        The format
        of the method's output correponds with the particle's algebra.

        Explain the output: # TODO
        Raises
        ------
         NotImplementedError
              Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def _initial_guess_builder_particle(self) -> None:
        """Constructs the initial guess for the differential equation
        solver that represents the particle phase behaviour.
        The format
        of the method's output correponds with the particle's algebra.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Particles Balance methods
    # ==================================================================
    @abstractmethod
    def get_kinetics(self) -> None:
        """Return kinetics from instance of reactor class

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def difussion_external_evaluation(self) -> None:
        """Method that evals external difusion mass coeficient
        and compares it with kinetics coeficient

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def catalyst_mass_balance(self) -> None:
        """Method that evals and returns the evaluated particle
        mass balances. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def catalyst_isothermal_energy_balance(self) -> None:
        """Method that evals and returns the evaluated particle energy balances
        in isothermal operation. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def catalyst_non_isothermal_energy_balance(self) -> None:
        """Method that evals and returns the evaluated particle energy balance
        in non-isothermal operation. The format of the method's returns
        correspond to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            "Abstract method not implemented."
        """
        raise NotImplementedError("Abstract method not implemented.")

    @abstractmethod
    def simulate(self) -> None:
        """Simulates the particle given the mass_balance_data,
        energy_balance_data, pressure_balance_data.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")
