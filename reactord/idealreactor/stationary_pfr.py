import numpy as np
from _collections_abc import Callable

from reactord import ReactorBase
from reactord.mix import AbstractMix
from reactord.utils import vectorize


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

        super().__init__(
            self,
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            options=options,
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

    # Mass settings
    def set_mass_balance_data(
        self,
        molar_flux_in: list[float],
        molar_flux_out: list[float],
        catalyst=None,
    ) -> None:
        """Method that recieves and instantiates the neccesary
        parameters to solve the mass balance in the reactor's bulk.

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
        catalyst : AbstractCatalyst, optional # TODO
           Stationary catalyst particle, by default None

        Raises
        ------
        IndexError
            molar_flux_in length must be equal to molar_flux_out
        IndexError
            molar_flux_in and molar_flux_out length must be equal to
            mixture's substance number.
        """

        # ==============================================================
        # Validation
        # ==============================================================

        if len(molar_flux_in) != len(molar_flux_out):
            raise IndexError(
                "molar_flux_in length must be equal to molar_flux_out"
            )

        if len(molar_flux_in) != len(self.mix):
            raise IndexError(
                "molar_flux_in and molar_flux_out length must be equal to "
                "mixture's substance number."
            )

        self._molar_flux_in = molar_flux_in
        self._molar_flux_out = molar_flux_out

        if catalyst is None:
            self._catalyst_operation = "homogeneous"
            self._mass_balance_func = self._homogeneous_mass_balance
        else:
            self._catalyst_operation = "heterogeneous"
            self._mass_balance_func = self._heterogeneous_mass_balance

    # Temperature settings
    def set_isothermic_operation(self, isothermic_temperature: float) -> None:
        """Sets isothermic operation data needed for energy balance.

        Parameters
        ----------
        isothermic_temperature : float
            Isothermic reactor's temperature. [K]
        """

        self._thermal_operation = "isothermic"
        self._isothermic_temperature = isothermic_temperature
        self._energy_balance_func = self._isothermic_energy_balance

    def set_adiabatic_operation(self):
        """Not implemented

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented.")

    def set_non_isothermic_operation(self):
        """Not implemented

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented")

    # Pressure settings
    def set_isobaric_operation(self, isobaric_pressure: float):
        """Sets isobaric operation data needed for energy balance.

        Parameters
        ----------
        isobaric_pressure : float
            Isobaric reactor's pressure. [Pa]
        """

        self._pressure_operation = "isobaric"
        self._isobaric_pressure = isobaric_pressure
        self._pressure_balance_func = self._isobaric_pressure_balance

    def set_non_isobaric_operation(self):
        """Not implemented.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented.")

    # ==================================================================
    # Solvers aditional data needed - general used methods
    # ==================================================================

    def _grid_builder(self, grid_size: int) -> np.ndarray[float]:
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

        Parameters
        ----------
        grid_size : int
            Size of the grid of the reactor's length independent variable.

        Returns
        -------
        np.ndarray[float]
            Grid for the grid of the reactor's length independent variable.
        """
        
        return np.linspace(
            self.reactor_dim_minmax[0], self.reactor_dim_minmax[1], grid_size
        )

    def _border_condition_builder(self, ya, yb) -> None:
        """Constructs the border conditions for the differential
        equations that represents the bulk phase behaviour. The format
        of the method's output correponds with the reactor's algebra.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """
        

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

    def _mass_balance(self, grid, molar_fluxes, temperature, pressure) -> None:
        """Method that evals and returns the evaluated reactor's bulk
        mass balances. The format of the method's returns corresponds
        to the specific solver needs.

        Explain the output: # TODO

        Raises
        ------
        NotImplementedError
            Abstract method not implemented.
        """

        if self._catalyst_operation == "":
            raise ValueError("set_mass_balance_data method first")
        else:
            return self._mass_balance_func()

    @vectorize(signature="(n),(n),(n)->()", excluded={0})
    def _energy_balance(
        self,
        grid: list[float],
        temperature: list[float],
        refrigerant_temperature: list[float],
    ) -> None:

        """_summary_

        Returns
        -------
        _type_
            _description_
        """

        return self._energy_balance_func(
            grid, temperature, refrigerant_temperature
        )

    def _pressure_balance(
        self,
    ) -> None:
        """method that resturns the derivative of pressure respecto to
        reactors length on each reactor's grid node.

        Parameters
        ----------
        z : list or numpy.ndarray[float]
            Reactor's length grid.
        temperature : list or numpy.ndarray[float]
            Temperature on each reactor's grid node. [K]
        pressure : list or numpy.ndarray[float]
            Pressure on each reactor's grid node. [Pa]

        Returns
        -------
        list[float]
            Derivative of pressure respect to reactors length on each
            reactor's grid node.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented yet.
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

    # ==================================================================
    # Specifics mass balances
    # ==================================================================

    def _homogeneous_mass_balance(self) -> None:
        """Not implemented.

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Abstract method not implemented.")

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

    def _isothermic_energy_balance(
        self,
        grid: list[float],
        temperature: list[float],
        temperature_refrigerant: list[float],
        pressure: list[float],
    ) -> np.ndarray[float]:
        """Returns de derivative of temperature respect to reactor's
        length for a isothermic PFR.

        # Latex math: TODO

        dT/dz = 0

        Parameters
        ----------
        grid : list or numpy.ndarray[float]
            Reactor's length grid.
        temperature : list or numpy.ndarray[float]
            Temperature on each reactor's grid node. [K]
        temperature_refrigerant : list[float]
            Refrigerant temperature on each reactor's grid node. [K]
        pressure : list or numpy.ndarray[float]
            Pressure on each reactor's grid node. [Pa]

        Returns
        -------
        np.ndarray[float]
            Derivative of temperature respect to reactors length on each
            reactor's grid node. [K/m]
        """
        return 0

    def _homogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _heterogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _homogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _heterogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    # ==================================================================
    # Specifics pressure balances
    # ==================================================================

    def _isobaric_pressure_balance(
        self,
        grid: list[float],
        temperature: list[float],
        pressure: list[float],
    ) -> np.ndarray[float]:
        """Not implemented yet.

        Parameters
        ----------
        grid : list or numpy.ndarray[float]
            Reactor's length grid.
        temperature : list or numpy.ndarray[float]
            Temperature on each reactor's grid node. [K]
        pressure : list or numpy.ndarray[float]
            Pressure on each reactor's grid node. [Pa]

        Not implemented yet.

        # TODO

        Returns
        -------
        np.ndarray[float]
            Derivative of pressure respect to reactors length on each
            reactor's grid node. [Pa/m]

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet.")

    def _non_isobaric_pressure_balance(
        self,
        grid: list[float],
        temperature: list[float],
        pressure: list[float],
    ) -> np.ndarray[float]:
        """Not implemented yet.

        Parameters
        ----------
        grid : list or numpy.ndarray[float]
            Reactor's length grid.
        temperature : list or numpy.ndarray[float]
            Temperature on each reactor's grid node. [K]
        pressure : list or numpy.ndarray[float]
            Pressure on each reactor's grid node. [Pa]

        Not implemented yet.

        # TODO

        Returns
        -------
        np.ndarray[float]
            Derivative of pressure respect to reactors length on each
            reactor's grid node. [Pa/m]

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet.")

    # ==================================================================
    # Specifics solvers
    # ==================================================================

    def _homogeneous_solver(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _heterogeneous_solver(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")
