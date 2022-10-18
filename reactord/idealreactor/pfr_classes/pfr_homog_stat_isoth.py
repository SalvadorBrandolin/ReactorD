import numpy as np
from scipy.integrate import solve_bvp

from reactord.kinetics import Kinetics
from reactord.mix import AbstractMix
from reactord.reactorbase import ReactorBase
from reactord.utils import vectorize


class PfrHomogStatIsoth(ReactorBase):
    def __init__(
        self,
        mix: AbstractMix,
        list_of_reactions: list,
        stoichiometry: list,
        reactor_dims_minmax: list,
        transversal_area: float,
        pressure: float,
        reactor_isothermic_temperature: float,
        reactor_f_in: list,
        reactor_f_out: list,
        kinetic_argument: str = "concentration",
        **options
    ):

        self.kinetic: Kinetics = Kinetics(
            list_of_reactions=list_of_reactions,
            mix=mix,
            stoichiometry=stoichiometry,
            kinetic_argument=kinetic_argument,
            **options
        )

        self.stoichiometry = np.array(stoichiometry)
        self.mix = mix
        self.list_of_reactions = list_of_reactions
        self.reactor_dims_minmax = reactor_dims_minmax
        self.transversal_area = transversal_area
        self.pressure = pressure
        self.temperature = reactor_isothermic_temperature
        self.reactor_f_in = np.array(reactor_f_in)
        self.reactor_f_out = np.array(reactor_f_out)

        # Index where border conditions are given
        self._in_index = np.invert(np.isnan(reactor_f_in))
        self._out_index = np.invert(np.isnan(reactor_f_out))

        self._in_index = np.argwhere(self._in_index).ravel()
        self._out_index = np.argwhere(self._out_index).ravel()

    def _grid_builder(self, grid_size: int):
        dim_array = np.linspace(
            self.reactor_dims_minmax[0], self.reactor_dims_minmax[1], grid_size
        )
        return dim_array

    def _border_condition_builder(self, ya, yb):
        bc = np.array([])

        for i in self._in_index:
            bc = np.append(bc, ya[i] - self.reactor_f_in[i])

        for j in self._out_index:
            bc = np.append(bc, yb[j] - self.reactor_f_out[i])

        return bc

    def _initial_guess_builder(self, grid_size):
        n_comp = self.kinetic.num_substances
        initial_guess = np.zeros([n_comp, grid_size])

        for i in self._in_index:
            initial_guess[i, :] = np.full(grid_size, self.reactor_f_in[i])

        for j in self._out_index:
            initial_guess[j, :] = np.full(grid_size, self.reactor_f_out[j])

        return initial_guess

    @vectorize(signature="()->()", excluded={0})
    def _mass_balance(self, substances_reaction_rates):
        dfi_dz = substances_reaction_rates * self.transversal_area
        return dfi_dz

    def _reactor_energy_balance(self, grid_size):
        pass

    def _pressure_balance(self, grid_size):
        pass

    def _particle_mass_balance(self):
        pass

    def _particle_energy_balance(self):
        pass

    def _refrigerant_energy_balance(self):
        pass

    def simulate(self, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0):
        def odesystem(z, variables):
            molar_fluxes = np.transpose(variables)

            ri_rates, reaction_rates = self.kinetic.kinetic_eval(
                molar_fluxes, self.temperature, self.pressure
            )

            df_dz = self._mass_balance(ri_rates).T

            return df_dz

        z = self._grid_builder(grid_size)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(grid_size)

        sol = solve_bvp(
            odesystem,
            bc,
            z,
            in_guess,
            tol=tol,
            max_nodes=max_nodes,
            verbose=verbose,
        )

        return sol
