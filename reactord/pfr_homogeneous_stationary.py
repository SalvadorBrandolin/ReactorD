from ReactorBase import ReactorBase
from kinetics import Kinetics
from Mix import Abstract_Mix
from scipy.integrate import solve_bvp
import numpy as np


class PFR_Homogeneous_Stationary(ReactorBase):
    
    def __init__(
        self, 
        mix: Abstract_Mix,
        list_of_reactions: list[function],
        stoichiometry: list[float],
        reactor_dims_minmax: list[float], 
        transversal_area: float,
        pressure: float, 
        reactor_t_operation: str, 
        reactor_f_in: list[float|str], 
        reactor_f_out: list[float|str],
        reactor_t_in: float|str, 
        reactor_t_out: float|str,
        refrigerant_t_operation: str=None,
        refrigerant_mix: Abstract_Mix=None,
        refrigerant_f_in: float=None,
        refrigerant_t_in: float|str=None, 
        refrigerant_t_out: float|str=None,
        refrigerant_pressure: float=None, 
        heat_exchange_disposition: str='cocurrent',       
        u: float=0,
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
        self.reactor_t_operation = reactor_t_operation
        self.refrigerant_t_operation = refrigerant_t_operation
        self.reactor_f_in = reactor_f_in
        self.reactor_f_out = reactor_f_out
        self.reactor_t_in = reactor_t_in
        self.reactor_t_out = reactor_t_out
        self.refrigerant_mix = refrigerant_mix
        self.refrigerant_f_in = refrigerant_f_in
        self.heat_exchange_disposition = heat_exchange_disposition
        self.refrigerant_t_in = refrigerant_t_in
        self.refrigerant_t_out = refrigerant_t_out
        self.refrigerant_pressure = refrigerant_pressure 
        self.u = u

        self._mass_balance = np.vectorize(
            self._mass_balance, 
            excluded='self', 
            signature='(n)->(n)')

        self._reactor_energy_balance = np.vectorize(
            self._reactor_energy_balance, 
            excluded='self', 
            signature='(),(n),(),(),(),(m)->()')

        self._refrigerant_energy_balance = np.vectorize(
            self._refrigerant_energy_balance,
            excluded='self',
            signature='(),(),()->()')
        
    def _grid_builder(
        self, 
        grid_size: int
    ):
        dim_array = np.linspace(self.reactor_dims_minmax[0], 
                                self.reactor_dims_minmax[1], 
                                grid_size)
        return dim_array

    def _border_condition_builder_deprecated(self, ya, yb):
        bc = np.array([])

        #Border conditions for reagents, products and inerts molar flux
        for i,(fin, fout) in enumerate(zip(self.reactor_f_in, 
                                           self.reactor_f_out)):
            if fin == 'var':
                bc = np.append(bc, yb[i] - fout)
            else:
                bc = np.append(bc, ya[i] - fin)
        
        #Border condition for the reactive mixture's temperature
        if self.reactor_t_operation == 'isothermal':
            bc = np.append(bc, ya[-3] - self.reactor_t_in)
        else:
            if self.reactor_t_in == 'var':
                bc = np.append(bc, yb[-3] - self.reactor_t_out)
            else:
                bc = np.append(bc, ya[-3] - self.reactor_t_in)
        
        #Border condition for pressure
        bc = np.append(bc, ya[-2] - self.pressure)
        
        #Border condition for refrigerant's temperature
        if self.refrigerant_t_operation == 'isothermal':
            bc = np.append(bc, ya[-1] - self.refrigerant_t_in)
        elif self.refrigerant_t_operation == 'non-isothermal':
            if self.refrigerant_t_in == 'var':
                bc = np.append(bc, yb[-1] - self.refrigerant_t_out)
            else:
                bc = np.append(bc, ya[-1] - self.refrigerant_t_in)
        else:
            bc = np.append(bc, ya[-1] - 0)

        return bc

    def _border_condition_builder(self, ya, yb): 
        in_conditions = np.append(
            (self.reactor_f_in, 
            (self.reactor_t_in, 
            self.pressure, 
            self.refrigerant_t_in)
        )







    def _initial_guess_builder(self, grid_size):
        n_comp = len(self.mix)
        in_guess = np.zeros([n_comp + 3, grid_size])
        
        #Guess for reagents, products and inerts molar flux
        for i,(fin, fout) in enumerate(zip(self.reactor_f_in, 
                                           self.reactor_f_out)):
            if fin == 'var':
                in_guess[i,:] = np.full(grid_size, fout)
            else:
                in_guess[i,:] = np.full(grid_size, fin)
        
        #Guess for the reactive mixture's temperature
        if self.reactor_t_operation == 'isothermal':
            in_guess[-3,:] = np.full(grid_size, self.reactor_t_in)
        else:
            if self.reactor_t_in == 'var':
                in_guess[-3,:] = np.full(grid_size, self.reactor_t_out)
            else:
                in_guess[-3,:] = np.full(grid_size, self.reactor_t_in)
        
        #Guess for reactor's pressure
        in_guess[-2,:] = np.full(grid_size, self.pressure)

        #Guess for refrigerant's temperature
        if self.refrigerant_t_operation == 'isothermal':
            in_guess[-1,:] = np.full(grid_size, self.refrigerant_t_in)
        elif self.refrigerant_t_operation == 'non-isothermal':
            if self.refrigerant_t_in == 'var':
                in_guess[-1,:] = np.full(grid_size, self.refrigerant_t_out)
            else:
                in_guess[-1,:] = np.full(grid_size, self.refrigerant_t_in)
        else:
            in_guess[-1,:] = np.full(grid_size, 0)

        return in_guess
                
    def _mass_balance(self, substances_reaction_rates):
        dfi_dz = substances_reaction_rates * self.transversal_area
        return dfi_dz
    
    def _reactor_energy_balance(
        self, grid_size, molar_fluxes, temperature, pressure,
        refrigerant_temperature, reaction_rates):
        
        if self.reactor_t_operation == 'isothermal':
            dt_dz = 0
            return dt_dz
        
        a = self.transversal_area
        u = self.u
        t = temperature
        ta = refrigerant_temperature
        p = pressure
        
        r_enthalpies = self.kinetics.reaction_enthalpies(t, p)
        mix_heat_capacity = self.mix.mix_heat_capacity(molar_fluxes, t, p)

        reactions_heat = np.dot(reaction_rates, r_enthalpies)
        total_molar_flux = np.sum(molar_fluxes, axis=0)
        
        if self.refrigerant_t_operation is None:
            dt_dz = (a
                     * (- reactions_heat)
                     / (total_molar_flux * mix_heat_capacity))
            return dt_dz
        else:
            dt_dz = (a 
                     * ((u * (ta - t)) - reactions_heat) 
                     / (total_molar_flux * mix_heat_capacity))
            return dt_dz
    
    def _pressure_balance(self, grid_size):
        dp_dz = np.full(grid_size, 0)
        return dp_dz

    def _refrigerant_energy_balance(
        self, grid_size, reactor_temperature, refrigerant_temperature
        ):
        
        if self.refrigerant_t_operation is None:
            dta_dz = 0
            return dta_dz
        elif self.reactor_t_operation == 'isothermal':
            dta_dz = 0
            return dta_dz
        else:
            refr_heat_capacity = self.refrigerant_mix.mix_heat_capacity(
                    self.refrigerant_f_in, ta, self.refrigerant_pressure
            )
            a = self.transversal_area
            u = self.u
            t = reactor_temperature
            ta = refrigerant_temperature
            f_ref = np.sum(self.refrigerant_f_in)

            dta_dz = a * u * (t - ta) / (f_ref *  refr_heat_capacity)
            return dta_dz

    def solve(self, grid_size=1000,  tol=0.001, max_nodes=1000, verbose=0):
        
        def odesystem(z,vars):
            fs = np.array(vars[0:len(vars)-3]).T
            t = np.array(vars[-3])
            p = np.array(vars[-2])
            ta = np.array(vars[-1])

            grid_solver = np.size(fs, axis=0)

            r_i, r_rates = self.kinetics.kinetics_eval(fs, t, p) 
            
            df_dz = self._mass_balance(r_i)
            
            dt_dz = self._reactor_energy_balance(
                grid_solver, fs, t, p, ta, r_rates)

            dp_dz = self._pressure_balance(grid_solver)

            dta_dz = self._refrigerant_energy_balance(grid_solver, t, ta)

            dvars_dz = np.vstack((df_dz.T, dt_dz, dp_dz, dta_dz))
            return dvars_dz

        z = self._grid_builder(grid_size)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(grid_size)
        pressure = np.full(grid_size, self.pressure, order='F')

        sol = solve_bvp(
            odesystem, bc, z, in_guess, tol=tol, 
            max_nodes=max_nodes,verbose=verbose
        )
        
        return sol