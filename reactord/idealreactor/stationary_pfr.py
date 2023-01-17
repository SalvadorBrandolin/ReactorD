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

    For instantiation is recommended to use the alternative class constructors:
    from_isothermic_isobaric
    from_isothermic_noisobaric
    from_adiabatic_isobaric
    from_adiabatic_noisobaric
    from_noisothermic_isobaric
    from_noisothermic_noisobaric

    For example:
    pfr = StationaryPFR.from_isothermic_isobaric(...)

    Parameters
    ----------
    mix : AbstractMix
        Mixture object.
    list_of_reactions : List[Callable]
        List of functions that evaluate the reaction rates, each one defined by
        the user with the following format:
        callable(concentration_unit: list[float], temperature: float
        ) -> float.
        Where concentration_unit refers to the units in which the arguments of
        the kinetic laws are expressed, for instance, concentrations or
        partial_pressures.
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
        'concentration': substance concentration. [mol/m³]
        'partial_pressure': substance partial pressure. [Pa]
    reactor_dim_minmax : List[float]
        List containing the minimum and maximum [min, max]
        boundaries of the reactor's length. E.g: a reactor modeled
        in the boundaries 0 m to 3 m has a
        reactor_dim_minmax = [0, 3]. [m]
    transversal_area : float
        Transversal area of the reactor. [m²]
    molar_flow_in : dict, optional
        Dictionary containing the known inlet molar flows of the
        substances, by default {}.  TODO
    molar_flow_out : dict, optional
        Dictionary containing the known outlet molar flows of the
        substances, by default {}. TODO
    catalyst_particle : _type_, optional
        CatalystParticle object, by default None TODO
    isothermic_temperature: float
        Temperature in isothermic operation [K]
    temperature_in_out: dict
        Temperature at inlet and outlet length in non isothermic or adiabatic
        operation. Eg: {'in': 500} or {'out': 800}. [K]
    refrigerant: AbstractMix
        Substances present in the refrigerant mix TODO
    refrigerant_molar_flow : float
        Refrigerant molar flow TODO
    refrigerant_temperature_in: float
        Inlet refrigerant temperature [K]
    refrigerant_constant_temperature: bool
        True is the refrigerant temperature do not change
        False is the refrigerant temperature change  TODO
    refrigerant_flow_arrangement: str
        Current or countercurrent
    exchanger_wall_material: Substance
        Substance with a conductivity constant defined
    correlation_heat_transfer: str
        Model to do the correlation heat transfer
    isobaric_pressure: float
        Pressure of isobaric operation
    pressure_in_out: dict
        Pressure at inlet and outlet length in non isothermic operation.
        Eg: {'in': 1013250} or {'out': 101325}. [Pa]
    pressure_loss_equation: str
        Type of equation to calculate the pressure loss. Options:
        "packed bed reactor" for Ergun's equation
        "gas phase reaction" for reactors without catalyst
    packed_bed_porosity: float
        Packed bed porosity float between 0 and 1.
    packed_bed_particle_diameter: float
        Diameter of particle of packed bed. [m]
    fanning_factor: float
        Friction factor for pipe without catalyst. This number depends on
        Reynolds and pipe roughness.
    """

    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        molar_flow_in: dict = None,
        molar_flow_out: dict = None,
        catalyst_particle: dict = None,
        isothermic_temperature: float = None,
        temperature_in_out: dict = None,
        refrigerant: AbstractMix = None,
        refrigerant_molar_flow: float = None,
        refrigerant_temperature_in: float = None,
        refrigerant_constant_temperature: bool = None,
        refrigerant_flow_arrangement: str = None,
        exchanger_wall_material: Substance = None,
        correlation_heat_transfer: str = None,
        isobaric_pressure: float = None,
        pressure_in_out: dict = None,
        pressure_loss_equation: str = None,
        packed_bed_porosity: float = None,
        packed_bed_particle_diameter: float = None,
        fanning_factor: float = None,
        reaction_enthalpies: List[float] = None,
    ) -> None:

        self._kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            reaction_enthalpies=reaction_enthalpies,
        )
        # =====================================================================
        # Specific reactor arguments
        # =====================================================================
        self._name = "StationaryPFR"
        self.reactor_dim_minmax = reactor_dim_minmax
        self.transversal_area = transversal_area

        # =====================================================================
        # Mass balance data
        # =====================================================================
        self.molar_flow_in = molar_flow_in
        self.molar_flow_out = molar_flow_out
        self.catalyst_particle = catalyst_particle

        # =====================================================================
        # Energy balance data
        # =====================================================================
        self.isothermic_temperature = isothermic_temperature
        self.temperature_in_out = temperature_in_out
        self.refrigerant = refrigerant
        self.refrigerant_molar_flow = refrigerant_molar_flow
        self.refrigerant_temperature_in = refrigerant_temperature_in
        self.refrigerant_constant_temperature = (
            refrigerant_constant_temperature
        )
        self.refrigerant_flow_arrangement = refrigerant_flow_arrangement
        self.exchanger_wall_material = exchanger_wall_material
        self.correlation_heat_transfer = correlation_heat_transfer

        # =====================================================================
        # Pressure balance data
        # =====================================================================
        self.isobaric_pressure = isobaric_pressure
        self.pressure_in_out = pressure_in_out
        self.pressure_loss_equation = pressure_loss_equation
        self.packed_bed_porosity = packed_bed_porosity
        self.packed_bed_particle_diameter = packed_bed_particle_diameter
        self.fanning_factor = fanning_factor

        # =====================================================================
        # Configure the reactor
        # =====================================================================
        self._set_catalyst_operation()
        self._set_thermal_operation()
        self._set_pressure_operation()

    # =========================================================================
    # Init constructors
    # =========================================================================

    @classmethod
    def from_isothermic_isobaric(
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
    ) -> "StationaryPFR":
        """Instantiate isothermic isobaric StationaryPFR.

        Parameters
        ----------
        mix : AbstractMix
            Mixture object.
        list_of_reactions : List[Callable]
            List of functions that evaluate the reaction rates, each one
            defined by the user with the following format:
            callable(concentration_unit: list[float], temperature: float
            ) -> float.
            Where concentration_unit refers to the units in which the arguments
            of the kinetic laws are expressed, for instance, concentrations or
            partial_pressures.
        stoichiometry : List[float]
            A matrix that represents the stoichiometry of the reactive system.
            Each row represents one reaction contained in the list_of_reactions
            parameter and each column represents a substance in the mix
            parameter. The stoichiometry matrix entrances are the
            stoichiometric coefficients of each substance in each reaction.
        kinetic_argument : str
            This argument is used to define how to evaluate the composition
            of the reactive mixture inside the reactor. Options:
            'concentration': substance concentration. [mol/m³]
            'partial_pressure': substance partial pressure. [Pa]
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max] boundaries of
            the reactor's length. E.g: a reactor modeled in the boundaries 0 m
            to 3 m has a reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Transversal area of the reactor. [m²]
        isothermic_temperature : float
            Reactor's temperature. [K]
        isobaric_pressure : float
            Reactor's pressure. [Pa]
        molar_flow_in : dict, optional
            Dictionary containing the known inlet molar flows of the
            substances, by default {}.  TODO
        molar_flow_out : dict, optional
            Dictionary containing the known outlet molar flows of the
            substances, by default {}. TODO
        catalyst_particle : _type_, optional
            CatalystParticle object, by default None TODO

        Returns
        -------
        StationaryPFR
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
    def from_isothermic_noisobaric(
        cls,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        isothermic_temperature: float,
        pressure_in_out: dict,
        pressure_loss_equation: str,
        packed_bed_porosity: float = None,
        packed_bed_particle_diameter: float = None,
        fanning_factor: float = None,
        molar_flow_in: dict = {},
        molar_flow_out: dict = {},
        catalyst_particle=None,
    ) -> "StationaryPFR":
        """Instantiate isothermic nonisobaric StationaryPFR.

        Parameters
        ----------
        mix : AbstractMix
            Mixture object.
        list_of_reactions : List[Callable]
            List of functions that evaluate the reaction rates, each one
            defined by the user with the following format:
            callable(concentration_unit: list[float], temperature: float
            ) -> float.
            Where concentration_unit refers to the units in which the arguments
            of the kinetic laws are expressed, for instance, concentrations or
            partial_pressures.
        stoichiometry : List[float]
            A matrix that represents the stoichiometry of the reactive system.
            Each row represents one reaction contained in the list_of_reactions
            parameter and each column represents a substance in the mix
            parameter. The stoichiometry matrix entrances are the
            stoichiometric coefficients of each substance in each reaction.
        kinetic_argument : str
            This argument is used to define how to evaluate the composition
            of the reactive mixture inside the reactor. Options:
            'concentration': substance concentration. [mol/m³]
            'partial_pressure': substance partial pressure. [Pa]
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max] boundaries of
            the reactor's length. E.g: a reactor modeled in the boundaries 0 m
            to 3 m has a reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Transversal area of the reactor. [m²]
        isothermic_temperature : float
            Reactor's temperature. [K]
        pressure_in_out: dict
            Pressure at inlet and outlet length in non isothermic operation.
            Eg: {'in': 1013250} or {'out': 101325}. [Pa]
        pressure_loss_equation: str
            Type of equation to calculate the pressure loss. Options:
            "packed bed reactor" for Ergun's equation
            "gas phase reaction" for reactors without catalyst
        packed_bed_porosity: float
            Packed bed porosity float between 0 and 1.
        packed_bed_particle_diameter: float
            Diameter of particle of packed bed. [m]
        fanning_factor: float
            Friction factor for pipe without catalyst. This number depends on
            Reynolds and pipe roughness.
        molar_flow_in : dict, optional
            Dictionary containing the known inlet molar flows of the
            substances, by default {}.  TODO
        molar_flow_out : dict, optional
            Dictionary containing the known outlet molar flows of the
            substances, by default {}. TODO
        catalyst_particle : _type_, optional
            CatalystParticle object, by default None TODO

        Returns
        -------
        StationaryPFR
            Instantiated isothermic noisobaric StationaryPFR.
        """
        isothermic_noisobaric_pfr = cls(
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
            pressure_in_out=pressure_in_out,
            pressure_loss_equation=pressure_loss_equation,
            packed_bed_porosity=packed_bed_porosity,
            packed_bed_particle_diameter=packed_bed_particle_diameter,
            fanning_factor=fanning_factor,
        )
        return isothermic_noisobaric_pfr

    @classmethod
    def from_adiabatic_isobaric(
        cls,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        temperature_in_out: dict,
        isobaric_pressure: float,
        molar_flow_in: dict = {},
        molar_flow_out: dict = {},
        catalyst_particle: dict = None,
        reaction_enthalpies: List[float] = None,
    ) -> "StationaryPFR":
        """Instantiate adiabatic isobaric StationaryPFR.

        Parameters
        ----------
        mix : AbstractMix
            Mixture object.
        list_of_reactions : List[Callable]
            List of functions that evaluate the reaction rates, each one
            defined by the user with the following format:
            callable(concentration_unit: list[float], temperature: float
            ) -> float.
            Where concentration_unit refers to the units in which the arguments
            of the kinetic laws are expressed, for instance, concentrations or
            partial_pressures.
        stoichiometry : List[float]
            A matrix that represents the stoichiometry of the reactive system.
            Each row represents one reaction contained in the list_of_reactions
            parameter and each column represents a substance in the mix
            parameter. The stoichiometry matrix entrances are the
            stoichiometric coefficients of each substance in each reaction.
        kinetic_argument : str
            This argument is used to define how to evaluate the composition
            of the reactive mixture inside the reactor. Options:
            'concentration': substance concentration. [mol/m³]
            'partial_pressure': substance partial pressure. [Pa]
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max] boundaries of
            the reactor's length. E.g: a reactor modeled in the boundaries 0 m
            to 3 m has a reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Transversal area of the reactor. [m²]
        temperature_in_out: dict
            Temperature at inlet and outlet length in non isothermic or
            adiabatic operation. Eg: {'in': 500} or {'out': 800}. [K]
        isobaric_pressure: float
            Pressure of isobaric operation. [Pa]
        molar_flow_in : dict, optional
            Dictionary containing the known inlet molar flows of the
            substances, by default {}.  TODO
        molar_flow_out : dict, optional
            Dictionary containing the known outlet molar flows of the
            substances, by default {}. TODO
        catalyst_particle : _type_, optional
            CatalystParticle object, by default None TODO
        """
        adiabatic_isobaric_pfr = cls(
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            reactor_dim_minmax=reactor_dim_minmax,
            transversal_area=transversal_area,
            molar_flow_in=molar_flow_in,
            molar_flow_out=molar_flow_out,
            catalyst_particle=catalyst_particle,
            temperature_in_out=temperature_in_out,
            isobaric_pressure=isobaric_pressure,
            reaction_enthalpies=reaction_enthalpies,
        )

        return adiabatic_isobaric_pfr

    @classmethod
    def from_adiabatic_noisobaric(
        cls,
        mix: AbstractMix,
        list_of_reactions: List[Callable],
        stoichiometry: List[float],
        kinetic_argument: str,
        reactor_dim_minmax: List[float],
        transversal_area: float,
        temperature_in_out: dict,
        pressure_in_out: dict,
        pressure_loss_equation: str,
        packed_bed_porosity: float = None,
        packed_bed_particle_diameter: float = None,
        fanning_factor: float = None,
        molar_flow_in: dict = {},
        molar_flow_out: dict = {},
        catalyst_particle=None,
    ) -> "StationaryPFR":
        """Instantiate adiabatic nonisobaric StationaryPFR.

        Parameters
        ----------
        mix : AbstractMix
            Mixture object.
        list_of_reactions : List[Callable]
            List of functions that evaluate the reaction rates, each one
            defined by the user with the following format:
            callable(concentration_unit: list[float], temperature: float
            ) -> float.
            Where concentration_unit refers to the units in which the arguments
            of the kinetic laws are expressed, for instance, concentrations or
            partial_pressures.
        stoichiometry : List[float]
            A matrix that represents the stoichiometry of the reactive system.
            Each row represents one reaction contained in the list_of_reactions
            parameter and each column represents a substance in the mix
            parameter. The stoichiometry matrix entrances are the
            stoichiometric coefficients of each substance in each reaction.
        kinetic_argument : str
            This argument is used to define how to evaluate the composition
            of the reactive mixture inside the reactor. Options:
            'concentration': substance concentration. [mol/m³]
            'partial_pressure': substance partial pressure. [Pa]
        reactor_dim_minmax : List[float]
            List containing the minimum and maximum [min, max] boundaries of
            the reactor's length. E.g: a reactor modeled in the boundaries 0 m
            to 3 m has a reactor_dim_minmax = [0, 3]. [m]
        transversal_area : float
            Transversal area of the reactor. [m²]
        temperature_in_out: dict
            Temperature at inlet and outlet length in non isothermic or
            adiabatic operation. Eg: {'in': 500} or {'out': 800}. [K]
        pressure_in_out: dict
            Pressure at inlet and outlet length in non isothermic operation.
            Eg: {'in': 1013250} or {'out': 101325}. [Pa]
        pressure_loss_equation: str
            Type of equation to calculate the pressure loss. Options:
            "packed bed reactor" for Ergun's equation
            "gas phase reaction" for reactors without catalyst
        packed_bed_porosity: float
            Packed bed porosity float between 0 and 1.
        packed_bed_particle_diameter: float
            Diameter of particle of packed bed. [m]
        fanning_factor: float
            Friction factor for pipe without catalyst. This number depends on
            Reynolds and pipe roughness.
        molar_flow_in : dict, optional
            Dictionary containing the known inlet molar flows of the
            substances, by default {}.  TODO
        molar_flow_out : dict, optional
            Dictionary containing the known outlet molar flows of the
            substances, by default {}. TODO
        catalyst_particle : _type_, optional
            CatalystParticle object, by default None TODO

        Returns
        -------
        StationaryPFR
            Instantiated adiabatic noisobaric StationaryPFR.
        """
        adiabatic_noisobaric_pfr = cls(
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            reactor_dim_minmax=reactor_dim_minmax,
            transversal_area=transversal_area,
            temperature_in_out=temperature_in_out,
            pressure_in_out=pressure_in_out,
            pressure_loss_equation=pressure_loss_equation,
            packed_bed_porosity=packed_bed_porosity,
            packed_bed_particle_diameter=packed_bed_particle_diameter,
            fanning_factor=fanning_factor,
            molar_flow_in=molar_flow_in,
            molar_flow_out=molar_flow_out,
            catalyst_particle=catalyst_particle,
        )

        return adiabatic_noisobaric_pfr

    @classmethod
    def from_noisothermic_isobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet")

    @classmethod
    def from_noisothermic_noisobaric(cls) -> None:
        """Not implemented yet.

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet")

    # =========================================================================
    # Abastract methods
    # Settings for mass, energy and pressure balances.
    # =========================================================================

    def _set_catalyst_operation(self) -> None:
        """Interpret the reactor's attributes to set the mass balance.

        Method that validates and interprets the reactor's attributes
        to select the correspondant solver and mass balance methods.
        The method is called in the reactor's instantiation.

        Raises
        ------
        ValueError
            Not border conditions specified for a substance molar flow.
        ValueError
            Two border conditions specified for a substance molar flow.
        """
        self._molar_flow_in_for_bc = []
        self._molar_flow_out_for_bc = []
        # Data validation
        #   All substances must have one and only one border condition:
        #   Flow data is stored in the two previous private arrays.
        #   The unexisting molar flux given are stored as None.
        #   The method _border_cond_and_initial_guesses handle with
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

        # Solver and mass balance setting
        if self.catalyst_particle is None:
            self._catalyst_operation = "homogeneous"
            self._mass_balance_func = self._homogeneous_mass_balance
            self._solver_func = self._homogeneous_solver
        else:
            self._catalyst_operation = "heterogeneous"
            self._mass_balance_func = self._heterogeneous_mass_balance
            self._solver_func = self._heterogeneous_solver

    def _set_thermal_operation(self) -> None:
        """Interpret the reactor's attributes to set the energy balance.

        Method that validates and interprets the reactor's attributes to select
        the correspondant solver and mass balance methods. The method is called
        in the reactor's instantiation.

        Raises
        ------
        ValueError
            Both isothermic and no isothermic data given becomes ambiguous.
        ValueError
            No border conditions given for temperature.
        NotImplementedError
            No isothermal operation not implemented yet.
        ValueError
            Ambiguous thermal specifications.
        ValueError
            No thermal specification were seted.
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
                    "because thermal operation becomes ambiguous."
                )
            else:
                self._thermal_operation = "isothermal"
                self._energy_balance_func = self._isothermic_energy_balance
                self._temperature_in_for_bc = self.isothermic_temperature
                self._temperature_out_for_bc = None
                self._refrigerant_temperature_in_for_bc = 0
                self._refrigerant_temperature_out_for_bc = None

        elif any(no_isothermic_args):
            if self.temperature_in_out is None:
                # Check if border conditions for temperature is given.
                raise ValueError(
                    "No border conditions given for temperature. "
                    "Specify temperature_in_out."
                )
            else:
                pass

            if self.temperature_in_out and not any(no_isothermic_args[1:]):
                # Check if temperature_in_out is given and no refrigerant or
                # heat exchange information are given. This defines and
                # adiabatic thermal operation.
                self._thermal_operation = "adiabatic"

                if self._catalyst_operation == "homogeneous":
                    self._energy_balance_func = (
                        self._homogeneous_adiabatic_energy_balance
                    )
                elif self._catalyst_operation == "heterogeneous":
                    self._energy_balance_func = (
                        self._heterogeneous_adiabatic_energy_balance
                    )

                self._temperature_in_for_bc = self.temperature_in_out.get("in")
                self._temperature_out_for_bc = self.temperature_in_out.get(
                    "out"
                )
                self._refrigerant_temperature_in_for_bc = 0
                self._refrigerant_temperature_out_for_bc = None

                # Init the standard reaction enthalpies
                self.kinetics.std_reaction_enthalpies_init()

            elif all(no_isothermic_args):
                raise NotImplementedError(
                    "No isothermal operation not implemented yet."
                )

            else:
                raise ValueError(
                    "Ambigous thermal specifications. \n"
                    "For an adiabatic operation specify only:\n"
                    "  - temperature_in_out\n"
                    "For a no isothermal operation, specify all:\n"
                    "  - temperature_in_out\n"
                    "  - refrigerant\n"
                    "  - refrigerant_molar_flow\n"
                    "  - refrigerant_temperature_in\n"
                    "  - refrigerant_constant_temperature\n"
                    "  - refrigerant_flow_arrangement\n"
                    "  - exchanger_wall_material\n"
                    "  - correlation_heat_transfer"
                )
        else:
            raise ValueError("No thermal specification were seted.")

    def _set_pressure_operation(self) -> None:
        """Interpret the reactor's attributes to set the energy balance.

        Method that validates and interprets the reactor's attributes
        to select the correspondant solver and mass balance methods.
        The method is called in the reactor's instantiation.

        Raises
        ------
        ValueError
            Both isobaric and non-isobaric configuration set is ambiguous.
        NotImplementedError
            Non-isobaric operation not implemented yet.
        ValueError
            No pressure configuration set.
        """
        no_isobaric_args = [
            self.pressure_in_out,
            self.pressure_loss_equation,
            self.packed_bed_porosity,
            self.packed_bed_particle_diameter,
            self.fanning_factor,
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
                    "    packed_bed_particle_diameter\n"
                    "    fanning_factor\n"
                    "because pressure operation becomes ambigous."
                )
            else:
                self._pressure_operation = "isobaric"
                self._pressure_balance_func = self._isobaric_pressure_balance
                self._pressure_in_for_bc = self.isobaric_pressure
                self._pressure_out_for_bc = None

        elif any(no_isobaric_args):
            # Check non-isobaric operation:
            if self.pressure_in_out is {}:
                raise ValueError(
                    "If noe isobaric operation is setting specify on border"
                    "condition for pressure"
                )
            else:
                self._pressure_in_for_bc = self.pressure_in_out.get("in")
                self._pressure_out_for_bc = self.pressure_in_out.get("out")
            if self.pressure_loss_equation == "packed bed reactor":
                self._pressure_balance_func = (
                    self._non_isobaric_pressure_balance_packed_bed_reactor
                )
                if (
                    self.packed_bed_porosity is None
                    or self.packed_bed_particle_diameter is None
                ):
                    raise ValueError(
                        "If non isobaric packed bed reactor is setting,"
                        "specify a packed bed porosity and particle diameter"
                    )
            elif self.pressure_loss_equation == "gas phase reaction":
                self._pressure_balance_func = (
                    self._non_isobaric_pressure_balance_gas_phase_reaction
                )
                if self.fanning_factor is None:
                    raise ValueError(
                        "If non isobaric gas phase reaction is setting,"
                        "specify a Fanning factor"
                    )
        else:
            raise ValueError("No pressure specification was set.")

    # =========================================================================
    # Solvers aditional data needed - general used methods
    # =========================================================================

    def _grid_builder(self, grid_size: int) -> np.ndarray:
        """Construct the reactor's length grid.

        Method to build the grid of the independent variables values.
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
            Size of the grid for the reactor's independent variable.

        Returns
        -------
        ndarray
            Grid for the reactor's independent variables.
        """
        return np.linspace(
            self.reactor_dim_minmax[0], self.reactor_dim_minmax[1], grid_size
        )

    def _border_cond_and_initial_guesses(
        self, grid_size: int
    ) -> Tuple[Callable, List[float]]:
        """Construct border condition and initial guess for solve_bvp.

        Constructs the border conditions for the differential
        equations that represents the bulk-phase behaviour and the
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
            Size of the grid for the reactor's independent variable.


        Returns
        -------
        tuple[Callable, List[float]]
            Border condition function needed and initial guess matrix for
            scipy.solve_bvp.
        """
        # ==============================================================
        # Border condition building
        # ==============================================================
        def border_conditions(ya: List[float], yb: List[float]) -> List[float]:
            """Border condition builder for scipy.solve_bvp.

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

        # =====================================================================
        # Initial guess building
        # =====================================================================

        initial_guess = np.zeros((len(self._inlet_information), grid_size))

        for idx, (inlet, outlet) in enumerate(
            zip(self._inlet_information, self._outlet_information)
        ):
            if inlet is None:
                initial_guess[idx, :] = np.full(grid_size, outlet)
            else:
                initial_guess[idx, :] = np.full(grid_size, inlet)

        return border_conditions, initial_guess

    # =========================================================================
    # Common reactor methods
    # =========================================================================
    def _refrigerant_energy_balance(
        self,
        length_coordinate: float,
        temperature: float,
        refrigerant_temperature: float,
    ) -> List[float]:
        """Evaluate the refrigerant energy balance.

        Explain more: TODO
        non_isothermic : TODO

        Parameters
        ----------
        length_coordinate : float
            Reactor's length coordinate.
        temperature : float
            Temperature at the length coordinate. [K]
        refrigerant_temperature : float
            Refrigerant's temperature at the length coordinate. [K]

        Returns
        -------
        float
            Derivative of the refrigerant's temperature respect to reactor's
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

    # =========================================================================
    # Specifics mass balances
    # =========================================================================

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
            Molar fluxes of each substance at the length_coordinate.
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

    # =========================================================================
    # Specific energy balances
    # =========================================================================

    def _isothermic_energy_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
        refrigerant_temperature: float,
    ) -> np.ndarray:
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
        pressure : float
            Pressure at the length coordinate. [Pa]
        refrigerant_temperature : float
            Refrigerant temperature at the length coordinate. [K]

        Returns
        -------
        ndarray
            Derivative of reactor's temperature respect to reactor's
            length. [K/m]
        """
        grid_length = np.size(length_coordinate)
        return np.zeros(grid_length)

    @vectorize(signature="(),(n),(),(),()->()", excluded={0})
    def _homogeneous_adiabatic_energy_balance(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
        refrigerant_temperature: float,
    ) -> np.ndarray:
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
        pressure : float
            Pressure at the length coordinate. [Pa]
        refrigerant_temperature : float
            Refrigerant temperature at the length coordinate. [K]

        Returns
        -------
        ndarray
            Derivative of reactor's temperature respect to reactor's
            length. [K/m]
        """
        a = self.transversal_area
        ri = self.kinetics.kinetic_eval(molar_fluxes, temperature, pressure)[1]
        dhi = self.kinetics.reaction_enthalpies(temperature, pressure)
        cp = self.mix.mix_heat_capacity(molar_fluxes, temperature, pressure)
        sum_fi = np.sum(molar_fluxes)

        dt_dz = a * np.dot(ri, -dhi) / (sum_fi * cp)

        return dt_dz

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

    # =========================================================================
    # Specific pressure balances
    # =========================================================================

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

    @vectorize(signature="(),(n),(),()->()", excluded={0})
    def _non_isobaric_pressure_balance_packed_bed_reactor(
        self,
        length_coordinate: float,
        molar_fluxes: List[float],
        temperature: float,
        pressure: float,
    ) -> np.ndarray:
        """Pressure balance evaluation for packed bed reactor operation.

        Pressure balance is calculated for each substance in mixture.

        Parameters
        ----------
        length_coordinate : float
            _description_
        molar_fluxes : List[float]
            _description_
        temperature : float
            _description_
        pressure : float
            _description_

        Returns
        -------
        np.ndarray
            _description_

        Raises
        ------
        NotImplementedError
            _description_
        """
        total_molar_flow = np.sum(molar_fluxes)
        mix_molecular_weight = self.mix.mixture_molecular_weight(molar_fluxes)
        total_mass_flow = total_molar_flow * mix_molecular_weight / 1000

        mass_density = self.mix.mass_density(
            molar_fluxes, temperature, pressure
        )

        mix_viscosity = self.mix.mixture_viscosity(
            molar_fluxes, temperature, pressure
        )

        f = total_mass_flow / self.transversal_area
        rho = mass_density
        phi = self.packed_bed_porosity
        dp = self.packed_bed_particle_diameter
        mu = mix_viscosity

        dp_dz = (
            -f
            / (rho * dp)
            * (1 - phi)
            / phi**3
            * (150 * (1 - phi) * mu / dp + 1.75 * f)
        )

        # import ipdb
        # ipdb.set_trace()
        return dp_dz

    def _non_isobaric_pressure_balance_gas_phase_reaction(self) -> List[float]:
        """Not implemented yet.

        TODO

        Raises
        ------
        NotImplementedError
            Not implemented yet.
        """
        raise NotImplementedError("Not implemented yet.")

    # =========================================================================
    # Specific solvers
    # =========================================================================
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

            mass_derivatives = self._mass_balance_func(
                length_coordinate, molar_fluxes, temperature, pressure
            ).T

            temperature_derivatives = self._energy_balance_func(
                length_coordinate,
                molar_fluxes,
                temperature,
                pressure,
                refrigerant_temperature,
            ).T

            pressure_derivatives = self._pressure_balance_func(
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
