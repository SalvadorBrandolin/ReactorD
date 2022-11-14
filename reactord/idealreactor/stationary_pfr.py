from _collections_abc import Callable

from reactord import Kinetics, ReactorBase
from reactord.mix import AbstractMix


class StationaryPFR(ReactorBase):
    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list[Callable],
        stoichiometry: list,
        kinetic_argument: str,
        reactor_dim_minmax: list[float],
        transversal_area: float,
        **options,
    ) -> None:

        self._kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            _not_reaction_enthalpies=True,
        )

        # ==============================================================
        # Specifics reactor arguments
        # ==============================================================

        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

    # ==================================================================
    # Abastract methods
    # Settings for mass, energy and pressure balance.
    # ==================================================================

    def set_mass_balance_data(
        self,
        molar_flux_in: list[float],
        molar_flux_out: list[float],
    ):
        """Method that recieves and instantiates the neccesary
        parameters to solve the mass balance in the reactor's bulk
        phase as attributes of the reactor.

        Parameters
        ----------
        molar_flux_in : list[float] or numpy.ndarray[float]
            List or or numpy.ndarray containing the known molar fluxes
            at the reactor's inlet. The ordering of the fluxes must be
            identical to the substance order in mixture. For unkown
            fluxes specify as numpy.nan. [mol/s]

        molar_flux_out : list[float] or numpy.ndarray[float]
            List or or numpy.ndarray containing the known molar fluxes
            at the reactor's outlet. The ordering of the fluxes must be
            identical to the substance order in mixture. For unkown
            fluxes specify as numpy.nan. [mol/s]
        """

        # ==============================================================
        # Validation
        # ==============================================================

        self._molar_flux_in = molar_flux_in
        self._molar_flux_out = molar_flux_out

    def set_energy_balance_data(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the energy balance in the reactor's bulk
        phase as attributes of the reactor. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    def set_pressure_balance_data(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the pressure balance in the reactor's bulk
        phase as attributes of the reactor. The method returns None.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Solvers aditional data needed - general used methods
    # ==================================================================

    def _grid_builder(self) -> None:
        """Method to build the grid of independent variables.
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

    def _border_condition_builder(self) -> None:
        """Constructs the border conditions for the differential
        equations that represents the bulk phase behaviour. The format
        of the method's output correponds with the reactor's algebra.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    def _initial_guess_builder(self) -> None:
        """Constructs the initial guess for the differential equation
        solver that represents the bulk phase behaviour.

        The format
        of the method's output correponds with the reactor's algebra.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

    # ==================================================================
    # Common reactors methods
    # ==================================================================

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

    def _reactor_energy_balance(self) -> None:
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
