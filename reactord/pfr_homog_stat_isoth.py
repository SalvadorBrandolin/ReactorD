from ReactorBase import ReactorBase
from kinetics import Kinetics
from Mix import Abstract_Mix
from scipy.integrate import solve_bvp
from decorators import vectorize
import numpy as np


class PFR_Homog_Stat_Isoth(ReactorBase):
    
    def __init__(
        self, 
        mix: Abstract_Mix,
        list_of_reactions: list[function],
        stoichiometry: list[float],
        reactor_dims_minmax: list[float], 
        transversal_area: float,
        pressure: float,
        reactor_isothermic_temperature: float, 
        reactor_f_in: list[float|str], 
        reactor_f_out: list[float|str],
        **options     
    ):

        ReactorBase.__init__(
            self,
            mix=mix,
            list_of_reactions=list_of_reactions,
            stoichiometry=stoichiometry,
            options=options 
        )

        self.reactor_dims_minmax = reactor_dims_minmax
        self.transversal_area = transversal_area
        self.pressure = pressure
        self.temperature = reactor_isothermic_temperature
        self.reactor_f_in = reactor_f_in
        self.reactor_f_out = reactor_f_out

        #Index where border conditions are given
        self._in_index = np.argwhere(
            np.isin(self.reactor_f_in, 'var', invert=True)).ravel()
        self._out_index = np.argwhere(
            np.isin(self.reactor_f_out,'var', invert=True)).ravel()

    def _grid_builder(self, grid_size: int):
        dim_array = np.linspace(self.reactor_dims_minmax[0], 
                                self.reactor_dims_minmax[1], 
                                grid_size)
        return dim_array

    def _border_condition_builder(self, ya, yb): 
        bc = np.array([])
        
        for i in self._in_index:
            bc = np.append(bc, ya[i] - self.reactor_f_in[i])

        for j in self._out_index:
            bc = np.append(bc, yb[j] - self.reactor_f_out[i])

        return bc

    def _initial_guess_builder(self, grid_size):
        n_comp = self.kinetic.total_substances
        initial_guess = np.zeros([n_comp, grid_size])
        
        for i in self._in_index:
            initial_guess[i, :] = np.full(grid_size, self.reactor_f_in[i])

        for j in self._out_index:
            initial_guess[j, :] = np.full(grid_size, self.reactor_f_out[j])

        return initial_guess

    @vectorize(signature='()->()', exclude='self')
    def _mass_balance(self, substances_reaction_rates):
        dfi_dz = substances_reaction_rates * self.transversal_area
        return dfi_dz
    
    def _reactor_energy_balance(self, grid_size):
        pass
    
    def _pressure_balance(self, grid_size):
        pass

    def _refrigerant_energy_balance(self):
        pass

    def solve(self, grid_size=1000,  tol=0.001, max_nodes=1000, verbose=0):
        
        def odesystem(z,vars):
            molar_fluxes = np.transpose(vars)

            ri_rates, reaction_rates = self.kinetic.kinetic_eval(
                molar_fluxes, self.temperature, self.pressure
            )

            df_dz = self._mass_balance(ri_rates)   

            return df_dz

        z = self._grid_builder(grid_size)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(grid_size)
        pressure = np.full(grid_size, self.pressure, order='F')

        sol = solve_bvp(
            odesystem, bc, z, in_guess, tol=tol, 
            max_nodes=max_nodes,verbose=verbose
        )
        
        return sol