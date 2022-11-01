from abc import ABCMeta, abstractmethod

from reactord.kinetics import Kinetics


class ReactorBase(metaclass=ABCMeta):
    """Abstract class interface for each reactor in ReactorD."""

    _kinetics: Kinetics = None
    _settings = {
        "reactor_type": "",
        "time_operation": "",
        "catalytic_operation": "",
        "thermal_operation": "",
        "pressure_operation": "",
    }

    # ==================================================================
    # Abastract methods
    # ODE/PDE reactors general used methods
    # ==================================================================

    @abstractmethod
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
        """
        raise NotImplementedError()

    @abstractmethod
    def _border_condition_builder(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def _initial_guess_builder(self) -> None:
        raise NotImplementedError()

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    @abstractmethod
    def _catalyst_mass_balance(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def _catalyst_energy_balance(self) -> None:
        raise NotImplementedError()

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    @property
    def kinetics(self) -> Kinetics:
        return self._kinetics

    @kinetics.setter
    def kinetics(self, new_kinetics: Kinetics) -> None:
        self._kinetics = new_kinetics

    @abstractmethod
    def _mass_balance(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def _reactor_energy_balance(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def _pressure_balance(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def _refrigerant_energy_balance(self) -> None:
        raise NotImplementedError()

    @abstractmethod
    def simulate(self) -> None:
        raise NotImplementedError()
