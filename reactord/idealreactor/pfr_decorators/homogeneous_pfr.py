from reactord.decoratorbase import DecoratorBase
from reactord.idealreactor.pfr_decorators.isothermic_pfr import IsothermicPFR
from reactord.reactorbase import ReactorBase
from reactord.utils import vectorize


class HomogeneousPFR(DecoratorBase):
    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor: ReactorBase = reactor

        if self._settings["time_operation"] == "stationary":
            self._massbalance_function = self._stationary_mass_balance
        else:
            self._massbalance_function = self._non_stationary_mass_balance

    # ==================================================================
    # Configuration methods: returns set_thermal(...(SpecificReactor))
    # ==================================================================

    def set_isothermic(self):

        self._settings["thermal_operation"] = "isothermic"

        return IsothermicPFR(self, isothermic_temperature)

    def set_adiabatic(self):

        self._settings["thermal_operation"] = "adiabatic"

        raise NotImplementedError("no implemented... yet")

    def set_non_isothermic(self):

        self._settings["thermal_operation"] = "non_isothermic"

        raise NotImplementedError("no implemented... yet")

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    def _initial_guess_builder(self, grid_size):
        ...

    # ==================================================================
    # Common reactors methods
    # ==================================================================

    @vectorize(signature="()->()", excluded={0})
    def _mass_balance(self, substances_reaction_rates):

        return self._massbalance_function(substances_reaction_rates)

    # ==================================================================
    # Specifics private methods
    # ==================================================================

    def _stationary_mass_balance(self, substances_reaction_rates):
        dfi_dz = substances_reaction_rates * self.transversal_area
        return dfi_dz

    def _non_stationary_mass_balance(self, substances_reaction_rates):
        # TODO
        raise NotImplementedError()
