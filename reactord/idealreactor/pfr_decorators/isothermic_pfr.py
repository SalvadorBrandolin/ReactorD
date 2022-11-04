from reactord.decoratorbase import DecoratorBase
from reactord.reactorbase import ReactorBase


class IsothermicPFR(DecoratorBase):
    def __init__(
        self, reactor: ReactorBase, isothermic_temperature: float
    ) -> None:

        self._reactor = reactor

        # Specifics attributes
        self.isothermic_temperature = isothermic_temperature

    # ==================================================================
    # Configuration methods: returns set_pressure(...(SpecificReactor))
    # ==================================================================

    def set_isobaric(self, isothermic_temperature: float):

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "isothermic",
                "pressure_operation": "isobaric",
            }
        )

        return IsobaricPFR(self, isothermic_temperature)

    def set_isobaric(self):

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "isothermic",
                "pressure_operation": "non-isobaric",
            }
        )

        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _grid_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._grid_builder(*args, **kargs)

    def _border_condition_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._border_condition_builder(
            *args, **kargs
        )

    def _initial_guess_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._initial_guess_builder(*args, **kargs)

    # ==================================================================
    # Heterogeneous reactors methods
    # ==================================================================

    def _catalyst_mass_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_mass_balance(*args, **kargs)

    def _catalyst_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_energy_balance(*args, **kargs)

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    def _mass_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._mass_balance(*args, **kargs)

    def _reactor_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._reactor_energy_balance(*args, **kargs)

    def _pressure_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._pressure_balance(*args, **kargs)

    def _refrigerant_energy_balance(self, *args, **kargs) -> None:
        return self._decorated_reactor._refrigerant_energy_balance(
            *args, **kargs
        )

    def simulate(self, *args, **kargs) -> None:
        return self._decorated_reactor.simulate(*args, **kargs)
