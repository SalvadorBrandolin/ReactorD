"""Stationary plug flow reactor.

Module containing the StationaryPFR class for stationary plug flow reactor.

Example
-------
TODO: example using the stationary plug flow reactor.
"""

from typing import Callable, List, Tuple

import numpy as np

from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase
from reactord.substance import Substance
from reactord.utils import vectorize

from scipy.integrate import solve_bvp


class StationaryPFR(ReactorBase):
    """Stationary plug flow reactor class.

    Parameters
    ----------
        mix : AbstractMix
            Mixture object.
        list_of_reactions : List[Callable]
            List containing functions to eval the reaction rates, each
            defined by the user with the form:
                callable(
                    concentration_unit: list[float], temperature: float
                ) -> float
            Where concentration_unit refers to the oncentration unit of
            measure that is argument of the kinetic law, eg:
            concentration [mol/m^3], partial pressure [Pa]. The return
            of the function must be single value with units: [mol/s/m^3]
        stoichiometry : List[float]
            Stoichiometry matrix of the reactive system. Each row
            represents each reaction contained in the list_of_reactions
            parameter, each column represents each substance in the mix
            parameter. The stoichiometry matrix entrances are the
            stoichiometric coefficients of each substance in each
            reaction. The ordering of the stochiometry matrix rows must
            be the same as the ordering of the reaction functions in
            list_of_reactions. Also, the order of the columns of the
            stoichiometry matrix must be the same as the substance
            ordering in the mix.
        kinetic_argument : str
            Kinetic argument used to eval the reaction defined by the
            user. Options:
                'concentration': substance concentration. [mol/m^3]
                'partial_pressure': substance partial pressure. [Pa]
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max]
            boundaries of the reactor's length. E.g: a reactor modeled
            in the boundaries 0 m to 3 m has a
            reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Tranversal area of the reactor. [m^2]
        molar_flow_in : List[float] or numpy.ndarray[float]
            List or or numpy.ndarray containing the known molar fluxes
            at the reactor's inlet. The ordering of the fluxes must be
            identical to the substance order in mixture. For unkown
            fluxes specify as numpy.nan. [mol/s]
        molar_flow_out : List[float] or numpy.ndarray[float]
            List or or numpy.ndarray containing the known molar fluxes
            at the reactor's outlet. The ordering of the fluxes must be
            identical to the substance order in mixture. For unkown
            fluxes specify as numpy.nan. [mol/s]
        catalyst : AbstractCatalyst, optional # TODO
           Stationary catalyst particle, by default None

    Attributes
    ----------
    reactor_dim_minmax : List[float]
        List containing the minimum and maximum [min, max] boundaries of
        the reactor's length. E.g: a reactor modeled in the boundaries 0
        m to 3 m has a reactor dim minmax = [0, 3]. [m]
    transversal_area : float
        Tranversal area of the reactor. [m^2]
    """

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        molar_flow_in: List[float],
        molar_flow_out: List[float],
        catalyst_particle=None,
        isothermic_temperature: float = None,
        temperature_in_out: dict = None,
        refrigerant: AbstractMix = None,
        refrigerant_molar_flow: float = None,
        refrigerant_temperature_in: float = None,
        refrigerant_constant_temperature: bool = False,
        refrigerant_flow_arrangement: str = "cocurrent",
        exchanger_wall_material: Substance = None,
        correlation_heat_transfer: str = "",
        isobaric_pressure: float = None,
        pressure_in_out: dict = None,
        pressure_loss_equation: str = "",
        packed_bed_porosity: float = 1,
        **options,
    ) -> None:

        super().__init__(
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            options=options,
        )

        # Mass balance arguments and settings:
        self.molar_flow_in = molar_flow_in
        self.molar_flow_out = molar_flow_out
        self.catalyst_particle = catalyst_particle
        
        self._set_catalyst_operation()

        # Specifics reactor arguments

        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

    # ==================================================================
    # Init constructors
    # ==================================================================

    @classmethod
    def set_isothermic_isobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    @classmethod
    def set_isothermic_noisobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    @classmethod
    def set_adiabatic_isobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    @classmethod
    def set_adiabatic_noisobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    @classmethod
    def set_noisothermic_isobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    @classmethod
    def set_noisothermic_noisobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotADirectoryError
            Not implemented yet.
        """
        raise NotADirectoryError("Not implemented yet")

    # ==================================================================
    # Abastract methods
    # Settings for mass, energy and pressure balance.
    # ==================================================================

    def _set_catalyst_operation(self) -> None:
        """_summary_

        Raises
        ------
        IndexError
            _description_
        IndexError
            _description_
        """
        if len(self.molar_flow_in) != len(self.molar_flow_out):
            raise IndexError(
                "molar_flow_in length must be equal to molar_flow_out"
            )

        if len(molar_flow_in) != len(self.mix):
            raise IndexError(
                "molar_flow_in and molar_flow_out length must be equal to "
                "mixture's substance number."
            )

        self._molar_flow_in = molar_flow_in
        self._molar_flow_out = molar_flow_out

        if catalyst_particle is None:
            self._catalyst_operation = "homogeneous"
            self._mass_balance_func = self._homogeneous_mass_balance
            self._solver_func = self._homogeneous_solver
        else:
            self._catalyst_operation = "heterogeneous"
            self._mass_balance_func = self._heterogeneous_mass_balance
            self._solver_func = self._heterogeneous_solver

    @abstractmethod
    def _set_thermal_operation(self):
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
    def _set_pressure_operation(self):
        """Method that recieves and instantiates the neccesary
        parameters to solve the isobaric pressure balance in the
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

    def _grid_builder(self, grid_size: int) -> List[float]:
        """Construct the reactor's length grid.

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

        Parameters
        ----------
        grid_size : int
            Size of the grid of the reactor's length independent variable.

        Returns
        -------
        List[float]
            Grid for the grid of the reactor's length independent variable.
        """
        return np.linspace(
            self.reactor_dim_minmax[0], self.reactor_dim_minmax[1], grid_size
        )

    def _border_cond_and_initial_guesses(
        self, grid_size: int
    ) -> Tuple[Callable, List[float]]:
        """Construct border condition and initial guess for solve_bvp.

        Constructs the border conditions for the differential
        equations that represents the bulk phase behaviour and the
        initial guess matrix for scipy.solve_bvp.

        Ordering of border conditions:
        [flux_1, flux_2, ... , flux_n, temp, pressure, refr_temp]

        temp: reactor's temperature
        refr_temp: refrigerant's temperature
        flux_i: molar flux of substance i of mix.
        pressure: reactor's pressure.

        Parameters
        ----------
        grid_size : int
            Size of the grid of the reactor's length independent variable.

        Returns
        -------
        tuple[Callable, List[float]]
            Border condition function needed and initial guess matrix for
            scipy.solve_bvp.

        Raises
        ------
        NotImplementedError
            Not implementd.
        NotImplementedError
            Not implementd.
        NotImplementedError
            Not implementd.
        NotImplementedError
            Not implementd.
        NotImplementedError
            Not implementd.
        ValueError
            Unkmown error. Maybe private attributes were changed?.
            Restart the reactor configuration.
        """
        # ==============================================================
        # Border condition building
        # ==============================================================
        def border_conditions(ya: List[float], yb: List[float]) -> List[float]:
            """Border condition for scipy.solve_bvp.

            Parameters
            ----------
            ya : List[float]
                Border conditions on reactor's inlet.
            yb : List[float]
                Border conditions on reactor's outlet.

            Returns
            -------
            List[float]
                Border conditions for scipy.solve_bvp.
            """
            bc = np.array([])

            for i in self._in_index:
                bc = np.append(bc, ya[i] - self._inlet_information[i])

            for j in self._out_index:
                bc = np.append(bc, yb[j] - self._outlet_information[j])

            return bc

        self._inlet_information = self._molar_flow_in
        self._outlet_information = self._molar_flow_out

        # Isothermic - isobaric

        if self._thermal_operation == "isothermic":
            if self._pressure_operation == "isobaric":

                # Asign the isothermic and isobaric conditions to inlet
                self._inlet_information = np.append(
                    self._inlet_information,
                    [
                        self._isothermic_temperature,
                        self._isobaric_pressure,
                        0.0,
                    ],
                )

                self._outlet_information = np.append(
                    self._outlet_information,
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                )

            elif self._pressure_operation == "non_isobaric":
                raise NotImplementedError()
                # TODO

        elif self._thermal_operation == "adiabatic":
            if self._pressure_operation == "isobaric":
                raise NotImplementedError()
                # TODO

            elif self._pressure_operation == "non_isobaric":
                raise NotImplementedError()
                # TODO

        elif self._thermal_operation == "non_isothermic":
            if self._pressure_operation == "isobaric":
                raise NotImplementedError()
                # TODO

            elif self._pressure_operation == "non_isobaric":
                raise NotImplementedError()
                # TODO

        else:
            raise ValueError("Thermal and pressure configurations failed.")

        self._in_index = np.invert(np.isnan(self._inlet_information))
        self._out_index = np.invert(np.isnan(self._outlet_information))

        self._in_index = np.argwhere(self._in_index).ravel()
        self._out_index = np.argwhere(self._out_index).ravel()

        # ==============================================================
        # Initial guess building
        # ==============================================================

        initial_guess = np.zeros((len(self._inlet_information), grid_size))

        for idx, (inlet, outlet) in enumerate(
            zip(self._inlet_information, self._outlet_information)
        ):
            if np.isnan(inlet):
                initial_guess[idx, :] = np.full(grid_size, outlet)
            else:
                initial_guess[idx, :] = np.full(grid_size, inlet)

        return border_conditions, initial_guess

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    @vectorize(signature="(),(n),(),()->(n)", excluded={0})
    def _mass_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
    ) -> List[float]:
        """Mass balance evaluation for each substance in mixture.

        latex math: # TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        List[float]
            Derivatives of mix substances' molar fluxes respect to
            reactor's length. [mol/s/m]

        Raises
        ------
        ValueError
            Triying to simulate without mass balance data setting.
            Configure mass balance data setting with the method
            set_mass_balance_data.
        """
        # TODO check the next if before, not on each mass balance call
        if self._catalyst_operation == "":
            raise ValueError("use set_mass_balance_data method first")
        else:
            return self._mass_balance_func(
                length_coordinate, molar_fluxes, temperature, pressure
            )

    @vectorize(signature="(),(n),(),(),()->()", excluded={0})
    def _energy_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        refrigerant_temperature: float,
        pressure,
    ) -> List[float]:
        """Energy balance evaluation for the reactor's bulk phase.

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        refrigerant_temperature : float
            Refrigerant temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        float
            Derivative of reactor's temperature respect to reactor's
            length. [K/m]

        Raises
        ------
        ValueError
            Triying to simulate without energy balance data setting.
            Configure energy balance data setting with the method
            set_isothermic_operation, set_adiabatic_operation or
            set_non_isothermic_operation.
        """
        # TODO check the next if before, not on each energy balance call
        if self._thermal_operation == "":
            raise ValueError(
                "use set_isothermic_operation, set_adiabatic_operation or"
                " set_non_isothermic_operation method first"
            )

        return self._energy_balance_func(
            length_coordinate,
            molar_fluxes,
            temperature,
            refrigerant_temperature,
            pressure,
        )

    @vectorize(signature="(),(n),(),()->()", excluded={0})
    def _pressure_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
    ) -> List[float]:
        """Pressure balance evaluation for each substance in mixture.

        latex math: # TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : list or List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        float
            Derivative of reactor's pressure respect to reactor's
            length. [mol/s/m]

        Raises
        ------
        ValueError
            Triying to simulate without pressure balance data setting.
            Configure pressurebalance data setting with the method
            set_isobaric_operation or set_non_isobaric_operation.
        """
        if self._pressure_operation == "":
            raise ValueError(
                "use set_isobaric_operation or set_non_isobaric_operation"
                " method first"
            )

        return self._pressure_balance_func(
            length_coordinate, molar_fluxes, temperature, pressure
        )

    @vectorize(signature="(),(),()->()", excluded={0})
    def _refrigerant_energy_balance(
        self,
        length_coordinate: float,
        temperature: float,
        refrigerant_temperature: float,
    ) -> List[float]:
        """Eval the refrigerant energy balance.

        Explain more: TODO
        non_isothermic : TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        temperature : float
            Temperature at the length coordinate. [K]
        refrigerant_temperature : float
            Refrigerant temperature at the length coordinate. [K]

        Returns
        -------
        float
            Derivative of refrigerant temperature respect to reactor's
            length.
        """
        if self._thermal_operation == "non_isothermic":
            raise NotImplementedError("Not implemented")
        else:
            return 0.0

    def simulate(
        self, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0
    ) -> List[float]:
        """Simulate the reactor.

        # TODO

        Parameters
        ----------
        grid_size : int, optional
            _description_, by default 1000
        tol : float, optional
            tolerance argument for scipy.solve_bvp , by default 0.001
        max_nodes : int, optional
            max_nodes argument for scipy.solve_bvp, by default 1000
        verbose : int, optional
            verbose argument for scipy.solve_bvp. Options 0 (False) or 1
            (True), by default 0

        Returns
        -------
        List[float]
            Solution matrix.
        """
        return self._solver_func(grid_size, tol, max_nodes, verbose)

    # ==================================================================
    # Specifics mass balances
    # ==================================================================

    def _homogeneous_mass_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
    ) -> List[float]:
        """Mass balance evaluation for each substance in mixture.

        latex math: # TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        List[float]
            Derivatives of mix substances' molar fluxes respect to
            reactor's length. [mol/s/m]
        """
        (
            self.substances_reaction_rates,
            self.reaction_rates,
        ) = self.kinetics.kinetic_eval(molar_fluxes, temperature, pressure)

        return self.substances_reaction_rates * self.transversal_area

    def _heterogeneous_mass_balance(self) -> None:
        """Not implemented yet.

        # TODO

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
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        refrigerant_temperature: float,
        pressure,
    ) -> List[float]:
        """Energy balance evaluation for the reactor's bulk phase.

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        refrigerant_temperature : float
            Refrigerant temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        float
            Derivative of reactor's temperature respect to reactor's
            length. [K/m]
        """
        return 0.0

    def _homogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented yet.

        TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _heterogeneous_adiabatic_energy_balance(self) -> None:
        """Not implemented yet.

        TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _homogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented yet.

        TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")

    def _heterogeneous_non_isothermic_energy_balance(self) -> None:
        """Not implemented yet.

        TODO

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
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
    ) -> List[float]:
        """Pressure balance evaluation for each substance in mixture.

        latex math: # TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        molar_fluxes : list or List[float]
            Molar flux of each substances at the length_coordinate.
            [mol/s]
        temperature : float
            Temperature at the length coordinate. [K]
        pressure : float
            Pressure at the length coordinate. [Pa]

        Returns
        -------
        float
            Derivative of reactor's pressure respect to reactor's
            length. [mol/s/m]
        """
        return 0.0

    def _non_isobaric_pressure_balance(self) -> List[float]:
        """Not implemented yet.

        TODO

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet.")

    # ==================================================================
    # Specifics solvers
    # ==================================================================

    def _homogeneous_solver(
        self, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0
    ) -> List[float]:
        """Simulate the homogeneous reactor.

        Parameters
        ----------
        grid_size : int, optional
            _description_, by default 1000
        tol : float, optional
            tolerance argument for scipy.solve_bvp , by default 0.001
        max_nodes : int, optional
            max_nodes argument for scipy.solve_bvp, by default 1000
        verbose : int, optional
            verbose argument for scipy.solve_bvp. Options 0 (False) or 1
            (True), by default 0

        Returns
        -------
        List[float]
            Solution matrix.
        """

        def odesystem(
            length_coordinate: List[float], variables: List[float]
        ) -> List[float]:
            """ODE system for scipy.solve_bvp.

            Parameters
            ----------
            length_coordinate : List[float]
                Reactor's length coordinate. [m]
            variables : List[float]
                Dependent variables on a two dimensional array in the
                order:
                    variables: [
                        [flux_1(z1), flux_1(z2), ..., flux_1(zn)],
                        [flux_2(z1), flux_2(z2), ..., flux_2(zn)],
                        .
                        .
                        .
                        [flux_m(z1), flux_m(z2), ..., flux_m(zn)],
                        [   T(z1)  ,   T(z2)   , ...,   T(zn)   ],
                        [   P(z1)  ,   P(z2)   , ...,   P(zn)   ],
                        [ Tref(z1) ,  Tref(z2) , ...,  Tref(zn) ],
                    ]
                where:
                    flux_i = molar flux of substance i
                    z_j = length_coordinate j
                    T = reactor's temperature
                    P = reactor's pressure
                    Tref = refrigerant's temperature

            Returns
            -------
            List[float]
                Derivatives matrix at each length_coordinate.
            """
            n_subs = len(self.mix)

            molar_fluxes = np.transpose(variables[0:n_subs, :])
            temperature = variables[-3, :]
            pressure = variables[-2, :]
            refrigerant_temperature = variables[-1, :]

            mass_derivatives = self._mass_balance(
                length_coordinate, molar_fluxes, temperature, pressure
            ).T

            temperature_derivatives = self._energy_balance(
                length_coordinate,
                molar_fluxes,
                temperature,
                refrigerant_temperature,
                pressure,
            ).T

            pressure_derivatives = self._pressure_balance(
                length_coordinate, molar_fluxes, temperature, pressure
            ).T

            refrigerant_temperature_derivatives = (
                self._refrigerant_energy_balance(
                    length_coordinate, temperature, refrigerant_temperature
                )
            ).T

            return np.vstack(
                (
                    mass_derivatives,
                    temperature_derivatives,
                    pressure_derivatives,
                    refrigerant_temperature_derivatives,
                )
            )

        grid = self._grid_builder(grid_size)

        bc, initial_guess = self._border_cond_and_initial_guesses(grid_size)

        simulation = solve_bvp(
            odesystem,
            bc,
            grid,
            initial_guess,
            tol=tol,
            max_nodes=max_nodes,
            verbose=verbose,
        )

        return simulation

    def _heterogeneous_solver(self) -> None:
        """Not implemented yet.

        # TODO

        Raises
        ------
        NotImplementedError
            Not implemented.
        """
        raise NotImplementedError("Not implemented yet.")
