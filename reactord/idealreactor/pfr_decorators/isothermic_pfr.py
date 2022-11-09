import numpy as np

from reactord.decoratorbase import DecoratorBase
from reactord.idealreactor.pfr_decorators.isobaric_pfr import IsobaricPFR
from reactord.reactorbase import ReactorBase


class IsothermicPFR(DecoratorBase):
    def __init__(
        self, reactor: ReactorBase, isothermic_temperature: float
    ) -> None:

        self._reactor = reactor

        # Specifics attributes
        self.isothermic_temperature = isothermic_temperature

        # ==============================================================
        # Reactor energy balance
        # ==============================================================

        if self._settings["time_operation"] == "stationary":
            if self._settings["catalytic_operation"] == "homogeneous":
                self._reactor_energy_balance_func = self.tuki1
            elif self._settings["catalytic_operation"] == "heterogeneous":
                self._reactor_energy_balance_func = self.tuki2

        elif self._settings["time_operation"] == "non_stationary":
            if self._settings["catalytic_operation"] == "homogeneous":
                self._reactor_energy_balance_func = self.tuki3
            elif self._settings["catalytic_operation"] == "heterogeneous":
                self._reactor_energy_balance_func = self.tuki4


        # ==============================================================
        # Catalyst energy balance
        # ==============================================================

        if self._settings["catalytic_operation"] == "homogeneous":
            pass
        else: 
            pass


    # ==================================================================
    # Configuration methods: returns set_pressure(...(SpecificReactor))
    # ==================================================================

    def set_isobaric(self, pressure: float):

        self._settings["pressure_operation"] = "isobaric"

        return IsobaricPFR(self, pressure)

    def set_non_isobaric(self):

        self._settings["pressure_operation"] = "non-isobaric"

        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _border_condition_builder(
        self, ya: np.ndarray, yb: np.ndarray
    ) -> np.ndarray:

        mass_bc = self._decorated_reactor._border_condition_builder(ya, yb)

        if self._settings["time_operation"] == "stationary":

            subs_num = len(self.kinetics.mix)

            reactor_temp_bc = ya[subs_num + 1] - self.isothermic_temperature
            refrigerant_temp_bc = ya[subs_num + 2] - 0

            temperature_bc = np.append(reactor_temp_bc, refrigerant_temp_bc)
            return np.append(mass_bc, temperature_bc)

        else:
            # TODO
            pass

    def _initial_guess_builder(self, grid_size) -> None:

        if self._settings["time_operation"] == "stationary":
            mass_guess = self._decorated_reactor._initial_guess_builder(
                grid_size
            )

        else:
            pass

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_energy_balance(self, *args, **kargs):
        if self._settings["catalytic_operation"] == "homogeneous":
            pass
        else: 
            self._settings["catalytic_operation"] == "heterogeneous"

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._reactor_energy_balance(*args, **kargs)

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._refrigerant_energy_balance(
            *args, **kargs
        )
