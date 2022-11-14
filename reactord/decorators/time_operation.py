from abc import ABCMeta, abstractmethod

from .decoratorbase import DecoratorBase


class TimeOperation(DecoratorBase, metaclass=ABCMeta):

    # ==================================================================
    # Configuration methods: returns set_catalysis(...(SpecificReactor))
    # ==================================================================

    @abstractmethod
    def set_homogeneous(self):
        raise NotImplementedError("abstractmethod not implemented")

    @abstractmethod
    def set_heterogeneous(self):
        raise NotImplementedError("abstractmethod not implemented")

    # ==================================================================
    # ReactorBase interface
    # ODE/PDE reactors general used methods
    # ==================================================================

    @abstractmethod
    def _grid_builder(self, *args, **kargs):
        """Specifics TimeOperation decorators must implement the
        _grid_builder method considering if the reactor is an ode
        or pde reactor. If algebraic reactor implement return None.

        Raises
        ------
        NotImplementedError
            abstractmethod not implemented. Specific decorator must
            implement the method.
        """

        raise NotImplementedError("abstractmethod not implemented")
