"""Stationary plug flow reactor.

Module containing the StationaryPFR class for stationary plug flow reactor.

Example
-------
TODO: example using the stationary plug flow reactor.
"""

from typing import Callable, List, Tuple

import numpy as np

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase
from reactord.substance import Substance
from reactord.utils import vectorize

from scipy.integrate import solve_bvp


class StationaryPFR(ReactorBase):
    """Stationary plug flow reactor class.

    For instantiation is recommended to use the setter class methods:
    set_isothermic_isobaric
    set_isothermic_noisobaric
    set_adiabatic_isobaric
    set_adiabatic_noisobaric
    set_noisothermic_isobaric
    set_noisothermic_noisobaric

    For example:
    pfr = StationaryPFR.set_isothermic_isobaric(...)

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
    reactor_dim_minmax : List[float]
        List containing the minimum and maximum [min, max]
        boundaries of the reactor's length. E.g: a reactor modeled
        in the boundaries 0 m to 3 m has a
        reactor_dim_minmax = [0, 3]. [m]
    transversal_area : float
        Transversal area of the reactor. [m^2]
    """

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        **kwargs,
    ) -> None:

        self._kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
        )
        # ==============================================================
        # Specifics reactor arguments
        # ==============================================================
        self._name = "StationaryPFR"
        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

        # ==============================================================
        # Mass balance data
        # ==============================================================
        self.molar_flow_in: dict = kwargs.get("molar_flow_in")
        self.molar_flow_out: dict = kwargs.get("molar_flow_out")
        self.catalyst_particle = kwargs.get("catalyst_particle")

        # ==============================================================
        # Energy balance data
        # ==============================================================
        self.isothermic_temperature: float = kwargs.get(
            "isothermic_temperature"
        )
        self.temperature_in_out: dict = kwargs.get("temperature_in_out")
        self.refrigerant: AbstractMix = kwargs.get("refrigerant")
        self.refrigerant_molar_flow: float = kwargs.get(
            "refrigerant_molar_flow"
        )
        self.refrigerant_temperature_in: float = kwargs.get(
            "refrigerant_temperature_in"
        )
        self.refrigerant_constant_temperature: bool = kwargs.get(
            "refrigerant_constant_temperature"
        )
        self.refrigerant_flow_arrangement: str = kwargs.get(
            "refrigerant_flow_arrangement"
        )
        self.exchanger_wall_material: Substance = kwargs.get(
            "exchanger_wall_material"
        )
        self.correlation_heat_transfer: str = kwargs.get(
            "correlation_heat_transfer"
        )

        # ==============================================================
        # Pressure balance data
        # ==============================================================
        self.isobaric_pressure: float = kwargs.get("isobaric_pressure")
        self.pressure_in_out: dict = kwargs.get("pressure_in_out")
        self.pressure_loss_equation: str = kwargs.get("pressure_loss_equation")
        self.packed_bed_porosity: float = kwargs.get("packed_bed_porosity")

        # ==============================================================
        # Configure the reactor
        # ==============================================================
        self._set_catalyst_operation()
        self._set_thermal_operation()
        self._set_pressure_operation()

    # ==================================================================
    # Init constructors
    # ==================================================================

    @classmethod
    def set_isothermic_isobaric(
        cls,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        isothermic_temperature: float,
        isobaric_pressure: float,
        molar_flow_in: dict = {},
        molar_flow_out: dict = {},
        catalyst_particle=None,
    ) -> ReactorBase:
        """Instantiate isothermic isobaric StationaryPFR.

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
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max]
            boundaries of the reactor's length. E.g: a reactor modeled
            in the boundaries 0 m to 3 m has a
            reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Transversal area of the reactor. [m^2]
        isothermic_temperature : float
            Reactor's temperature. [K]
        isobaric_pressure : float
            Reactor's pressure. [Pa]
        molar_flow_in : dict, optional
            Dictionary containing the known inlet molar flow of the
            substances, by default {}.  TODO
        molar_flow_out : dict, optional
            Dictionary containing the known outlet molar flow of the
            substances, by default {}. TODO
        catalyst_particle : _type_, optional
            CatalystParticle object, by default None TODO

        Returns
        -------
        ReactorBase
            Instantiated isothermic isobaric StationaryPFR.
        """
        isothermic_isobaric_pfr = cls(
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            reactor_dim_minmax=reactor_dim_minmax,
            transversal_area=transversal_area,
            molar_flow_in=molar_flow_in,
            molar_flow_out=molar_flow_out,
            catalyst_particle=catalyst_particle,
            isothermic_temperature=isothermic_temperature,
            isobaric_pressure=isobaric_pressure,
        )
        return isothermic_isobaric_pfr

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
        """Interpret the reactor's attributes to set the mass balance.

        Method that validates and interprets the reactor's attributes
        to select the correspondant solver and mass balance methods.
        The method is called in the reactor's instantiation.

        Raises
        ------
        ValueError
            Not border condition specified for a subtance molar flow.
        ValueError
            Two border condition specified for a subtance molar flow.
        """
        self._molar_flow_in_for_bc = []
        self._molar_flow_out_for_bc = []
        # Data validation
        #   All substances must have one and only one border condition:
        #   Flow data is stored in the two previous private arrays.
        #   The unexisting molar flux given are stored as None.
        #   the method _border_cond_and_initial_guesses handle with
        #   those two arrays.
        for substance in self.mix.substances:
            flow_inlet = self.molar_flow_in.get(substance.name)
            flow_outlet = self.molar_flow_out.get(substance.name)

            if (flow_inlet is None) and (flow_outlet is None):
                raise ValueError(
                    "Not molar flux (in or out) specified for substance: "
                    f"{substance.name}"
                )

            if (flow_inlet is not None) and (flow_outlet is not None):
                raise ValueError(
                    "Specify only molar_flow_in or molar_flow_out for "
                    f"substance: {substance.name}"
                )

            self._molar_flow_in_for_bc.append(flow_inlet)
            self._molar_flow_out_for_bc.append(flow_outlet)

        # Solver and mas balance setting
        if self.catalyst_particle is None:
            self._catalyst_operation = "homogeneous"
            self._mass_balance_func = self._homogeneous_mass_balance
            self._solver_func = self._homogeneous_solver
        else:
            self._catalyst_operation = "heterogeneous"
            self._mass_balance_func = self._heterogeneous_mass_balance
            self._solver_func = self._heterogeneous_solver

    def _set_thermal_operation(self):
        """Interpret the reactor's attributes to set the energy balance.

        Method that validates and interprets the reactor's attributes
        to select the correspondant solver and mass balance methods.
        The method is called in the reactor's instantiation.

        Raises
        ------
        ValueError
            Both isothermic and no isothermic configuration set is
            ambiguous.
        NotImplementedError
            No isothermic operation implemented yet.
        ValueError
            No thermal configuration seted.
        """
        no_isothermic_args = [
            self.temperature_in_out,
            self.refrigerant,
            self.refrigerant_molar_flow,
            self.refrigerant_temperature_in,
            self.refrigerant_constant_temperature,
            self.refrigerant_flow_arrangement,
            self.exchanger_wall_material,
            self.correlation_heat_transfer,
        ]

        if self.isothermic_temperature is not None:
            # Check isothermal operation:
            if any(no_isothermic_args):
                raise ValueError(
                    "If isothermic_temperature is specified, don't specify"
                    " any:\n"
                    "    temperature_in_out\n"
                    "    refrigerant\n"
                    "    refrigerant_molar_flow\n"
                    "    refrigerant_temperature_in\n"
                    "    refrigerant_constant_temperature\n"
                    "    refrigerant_flow_arrangement\n"
                    "    exchanger_wall_material\n"
                    "    correlation_heat_transfer\n"
                    "because thermal operation becomes ambigous."
                )
            else:
                self._thermal_operation = "isothermal"
                self._energy_balance_func = self._isothermic_energy_balance
                self._temperature_in_for_bc = self.isothermic_temperature
                self._temperature_out_for_bc = None
                self._refrigerant_temperature_in_for_bc = 0
                self._refrigerant_temperature_out_for_bc = None

        elif any(no_isothermic_args):
            # Check non isothermal operation:
            raise NotImplementedError(
                "No isothermic operation not implemented yet"
            )
        else:
            raise ValueError("No thermal specification where seted.")

    def _set_pressure_operation(self):
        """Interpret the reactor's attributes to set the energy balance.

        Method that validates and interprets the reactor's attributes
        to select the correspondant solver and mass balance methods.
        The method is called in the reactor's instantiation.

        Raises
        ------
        ValueError
            Both isobaric and no isobaric configuration set is
            ambiguous.
        NotImplementedError
            No isobaric operation implemented yet.
        ValueError
            No pressure configuration seted.
        """
        no_isobaric_args = [
            self.pressure_in_out,
            self.pressure_loss_equation,
            self.packed_bed_porosity,
        ]

        if self.isobaric_pressure is not None:
            # Check isobaric operation:
            if any(no_isobaric_args):
                raise ValueError(
                    "If isobaric_pressure is specified, don't specify"
                    " any:\n"
                    "    self.pressure_in_out\n"
                    "    self.pressure_loss_equation\n"
                    "    packed_bed_porosity\n"
                    "because pressure operation becomes ambigous."
                )
            else:
                self._pressure_operation = "isobaric"
                self._pressure_balance_func = self._isobaric_pressure_balance
                self._pressure_in_for_bc = self.isobaric_pressure
                self._pressure_out_for_bc = None

        elif any(no_isobaric_args):
            # Check no isobaric operation:
            raise NotImplementedError(
                "No isobaric operation not implemented yet"
            )
        else:
            raise ValueError("No pressure specification where seted.")

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

        self._inlet_information = np.append(
            self._molar_flow_in_for_bc,
            (
                self._temperature_in_for_bc,
                self._pressure_in_for_bc,
                self._refrigerant_temperature_in_for_bc,
            ),
        )
        self._outlet_information = np.append(
            self._molar_flow_out_for_bc,
            (
                self._temperature_out_for_bc,
                self._pressure_out_for_bc,
                self._refrigerant_temperature_out_for_bc,
            ),
        )

        self._in_index = np.argwhere(
            np.not_equal(self._inlet_information, None)
        ).ravel()
        self._out_index = np.argwhere(
            np.not_equal(self._outlet_information, None)
        ).ravel()

        # ==============================================================
        # Initial guess building
        # ==============================================================

        initial_guess = np.zeros((len(self._inlet_information), grid_size))

        for idx, (inlet, outlet) in enumerate(
            zip(self._inlet_information, self._outlet_information)
        ):
            if inlet is None:
                initial_guess[idx, :] = np.full(grid_size, outlet)
            else:
                initial_guess[idx, :] = np.full(grid_size, inlet)

        return border_conditions, initial_guess

    # ==================================================================
    # Common reactors methods
    # ==================================================================

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
        return self._mass_balance_func(
            length_coordinate, molar_fluxes, temperature, pressure
        )

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
        grid_length = np.size(length_coordinate)
        return np.zeros(grid_length)

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

    @vectorize(signature="(),(n),(),()->(n)", excluded={0})
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
        grid_length = np.size(length_coordinate)
        return np.zeros(grid_length)

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
        grid_length = np.size(length_coordinate)
        return np.zeros(grid_length)

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
