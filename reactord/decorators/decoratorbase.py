from dataclasses import dataclass

from reactord.reactorbase import ReactorBase


@dataclass
class DecoratorBase(ReactorBase):
    """
    Decorators interface
    """

    _reactor: ReactorBase = None
    _initialized: bool = False
    #     Internal variable used in the __setattr__ method. Indicates if
    #     the object was alredy instantiated or not. The value changes
    #     true in the __post_init method.

    def __init__(self, reactor: ReactorBase) -> None:
        self._reactor: ReactorBase = reactor

    def __post_init__(self):
        """Common __post_init__ method for all decorator classes."""
        self._initialized = True

    @property
    def _decorated_reactor(self) -> ReactorBase:
        return self._reactor

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

    def _catalyst_border_conditions(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_border_conditions(
            *args, **kargs
        )

    def _catalyst_initial_guess_builder(self, *args, **kargs) -> None:
        return self._decorated_reactor._catalyst_initial_guess_builder(
            *args, **kargs
        )

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

    # ==================================================================
    # Dunders
    # ==================================================================

    def __getattr__(self, name: str):
        """Custom __getattr__ method for all decorators. If the asked
        attributed is not found in the actual object, then is searched
        in the decorated object.

        Parameters
        ----------
        name : str
            name of the attribute

        Returns
        -------
        Asked Attribute type
            Attribute value
        """
        return getattr(self._decorated_reactor, name)

    def __setattr__(self, name: str, value):
        """Custom __setattr__ methods for all decorators. Uses the
        private variable _initialized to know if the object was
        instantiated or not. If not, the default __setattr__ is used.
        If the object is already instantiated a custom __setattr__ is
        called:
            If the attribute to set is not in the object, the
            __setattr__ method of the decorated reactor is called.
            Else the default __setattr__ is used.

        Parameters
        ----------
        name : str
            Name of the attribute
        value : Any
            Attribute asignation value.
        """
        if self._initialized:
            if not (name in self.__dict__):
                self._decorated_reactor.__setattr__(name, value)
            else:
                object.__setattr__(self, name, value)
        else:
            object.__setattr__(self, name, value)
