from abc import ABCMeta, abstractmethod


class ReactorBase(metaclass=ABCMeta):
    @abstractmethod
    def _grid_builder(self):
        pass

    @abstractmethod
    def _border_condition_builder(self):
        pass

    @abstractmethod
    def _initial_guess_builder(self):
        pass

    @abstractmethod
    def _mass_balance(self):
        pass

    @abstractmethod
    def _reactor_energy_balance(self):
        pass

    @abstractmethod
    def _pressure_balance(self):
        pass

    @abstractmethod
    def _particle_mass_balance(self):
        pass

    @abstractmethod
    def _particle_energy_balance(self):
        pass

    @abstractmethod
    def _refrigerant_energy_balance(self):
        pass

    @abstractmethod
    def simulate(self):
        pass
