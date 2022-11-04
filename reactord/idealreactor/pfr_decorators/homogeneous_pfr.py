from reactord.decoratorbase import DecoratorBase
from reactord.idealreactor.pfr_decorators.isothermic_pfr import IsothermicPFR
from reactord.reactorbase import ReactorBase
from reactord.utils import vectorize


class HomogeneousPFR(DecoratorBase):
    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor: ReactorBase = reactor

    # ==================================================================
    # Configuration methods: returns set_thermal(...(SpecificReactor))
    # ==================================================================

    def set_isothermic(self, isothermic_temperature: float):

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "isothermic",
                "pressure_operation": "",
            }
        )

        return IsothermicPFR(self, isothermic_temperature)

    def set_adiabatic(self):

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "isothermic",
                "pressure_operation": "",
            }
        )

        raise NotImplementedError("no implemented... yet")

    def set_non_isothermic(self):

        self._settings = dict(
            {
                "reactor_type": "Piston flow reactor (PFR)",
                "time_operation": "stationary",
                "catalytic_operation": "homogeneous",
                "thermal_operation": "isothermic",
                "pressure_operation": "",
            }
        )

        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    @vectorize(signature="()->()", excluded={0})
    def _mass_balance(self, substances_reaction_rates):
        dfi_dz = substances_reaction_rates * self.transversal_area
        return dfi_dz
