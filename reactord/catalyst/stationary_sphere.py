import numpy as np


from abstract_catalyst import AbstractCatalyst
from reactord.idealreactor.stationary_pfr import StationaryPFR
from reactord.utils import vectorize


class StationarySphere(AbstractCatalyst):
    def __init__(
        self,
        stationary_pfr: StationaryPFR,
        particle_radio: float,
        particle_density: float = None,
        porosity: float = None,
        total_particle_mass: float = None,
        superficial_area: float = None,
        particle_superficial_concentration: list = None,
        reactor_isothermic_temperature: float = None,
        particle_superficial_temperature: float = None,
    ):

        self.stationary_pfr = stationary_pfr
        self.particle_radio = particle_radio
        self.particle_density = particle_density
        self.porosity = porosity
        self.total_particle_mass = total_particle_mass
        self.superficial_area = superficial_area
        self.particle_superficial_concentration = (
            particle_superficial_concentration
        )
        self.reactor_isothermic_temperature = reactor_isothermic_temperature
        self.particle_superficial_temperature = (
            particle_superficial_temperature
        )

    def _grid_builder_particle(self, grid_size: int):
        dim_array = np.linspace(0, self.particle_radio, grid_size)
        return dim_array

    def _border_condition_builder_particle(self, ya, yb):
        bc = np.array([])

        for i in self._in_index:
            bc = np.append(
                bc, ya[i] - self.particle_superficial_concentration[i]
            )

        return bc

    def _initial_guess_builder_particle(self, grid_size):
        n_comp = self.kinetic.num_substances
        initial_guess = np.zeros([n_comp, grid_size])

        for i in self._in_index:
            initial_guess[i, :] = np.full(
                grid_size, self.particle_superficial_concentration[i]
            )

        return initial_guess

    def _get_kinetics(self):
        pass

    @vectorize(signature="()->()", excluded={0})
    def _catalyst_mass_balance(self):
        pass

    def _catalyst_energy_balance(self):
        pass

    def simulate(self, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0):
        pass
