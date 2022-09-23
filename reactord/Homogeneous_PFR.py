#from ReactorBase import ReactorBase
from scipy.integrate import solve_bvp
import numpy as np


class Homogeneous_PFR:
    
    def __init__(
        self, mix, kinetic, reactor_dims_minmax, transversal_area,
        pressure, reactor_t_operation, 
        reactor_f_in, reactor_f_out,
        reactor_t_in, reactor_t_out,
        refrigerant_t_operation=None,
        refrigerant_mix=None,
        refrigerant_f_in=None,
        refrigerant_t_in=None, refrigerant_t_out=None, 
        heat_exchange_disposition='cocurrent',       
        u=0        
    ):

        self.mix = mix
        self.kinetic = kinetic
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
        self.u = u
        
    def _grid_builder(self, grid_size):
        dim_array = np.linspace(self.reactor_dims_minmax[0], 
                                self.reactor_dims_minmax[1], 
                                grid_size)
        return dim_array

    def _border_condition_builder(self, ya, yb):
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
            bc = np.append(bc, ya[-2] - self.reactor_t_in)
        else:
            if self.reactor_t_in == 'var':
                bc = np.append(bc, yb[-2] - self.reactor_t_out)
            else:
                bc = np.append(bc, ya[-2] - self.reactor_t_in)
        
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

    def _initial_guess_builder(self, grid_size):
        n_comp = len(self.mix)
        in_guess = np.zeros([n_comp + 2, grid_size])
        
        #Guess for reagents, products and inerts molar flux
        for i,(fin, fout) in enumerate(zip(self.reactor_f_in, 
                                           self.reactor_f_out)):
            if fin == 'var':
                in_guess[i,:] = np.full(grid_size, fout)
            else:
                in_guess[i,:] = np.full(grid_size, fin)
        
        #Guess for the reactive mixture's temperature
        if self.reactor_t_operation == 'isothermal':
            in_guess[-2,:] = np.full(grid_size, self.reactor_t_in)
        else:
            if self.reactor_t_in == 'var':
                in_guess[-2,:] = np.full(grid_size, self.reactor_t_out)
            else:
                in_guess[-2,:] = np.full(grid_size, self.reactor_t_in)
        
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
        self, grid_size, molar_fluxes, reactor_temperatures, mix_heat_capacity,
        refrigerant_temperatures, reaction_rates, reaction_enthalpies):
        
        if self.reactor_t_operation == 'isothermal':
            dt_dz = np.zeros(grid_size)
            return dt_dz
        
        a = self.transversal_area
        u = self.u
        t = reactor_temperatures
        ta = refrigerant_temperatures
        reactions_heat = np.dot(reaction_rates, reaction_enthalpies)
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
    
    def _refrigerant_energy_balance(
        self, grid_size, reactor_temperature, refrigerant_temperature, 
        refrigerant_heat_capacity
        ):
        
        if self.refrigerant_t_operation is None:
            dta_dz = np.zeros(grid_size)
            return dta_dz
        elif self.reactor_t_operation == 'isothermal':
            dta_dz = np.zeros(grid_size)
            return dta_dz
        else:
            a = self.transversal_area
            u = self.u
            t = reactor_temperature
            ta = refrigerant_temperature
            f_ref = np.sum(self.refrigerant_f_in)

            dta_dz = a * u * (t - ta) / (f_ref * refrigerant_heat_capacity)
            return dta_dz

    def solve(self, grid_size=1000,  tol=0.001, max_nodes=1000, verbose=0):
        
        def odesystem(z,vars):
            fs = np.array(vars[0:len(vars)-1])
            t_reactor = np.array(vars[-2])
            t_refrigerant = np.array(vars[-1])
            
            r_i, r_rates = kinetic_eval(fs, t_reactor, pressure) 
            
            r_enthalpies_eval = reaction_enthalpies(t_reactor, pressure)

            reactor_heat_capacities = mix_heat_capacity(fs, t_reactor)

            refrig_heat_capacities = refr_heat_capacities(
                t_refrigerant, pressure
            )
            
            df_dz = self._mass_balance(r_i)
            
            dt_dz = self._reactor_energy_balance(
                grid_size, fs, t_reactor, reactor_heat_capacities,
                t_refrigerant, r_rates, r_enthalpies_eval
            )

            dta_dz = self._refrigerant_energy_balance(
                grid_size, t_reactor, t_refrigerant, refrig_heat_capacities
            )

            
            dvars_dz = np.vstack(df_dz, dt_dz, dta_dz)
            return dvars_dz

        z = self._grid_builder(grid_size)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(grid_size)
        pressure = np.full(grid_size, self.pressure)

        #Function vectorization
        kinetic_eval = np.vectorize(
            self.kinetic.kinetic_eval, signature='(n),(m),(m)->(m),(m)')
       
        reaction_enthalpies = np.vectorize(self.kinetic.reaction_enthalpies)
        
        mix_heat_capacity = np.vectorize(
            self.mix.mix_heat_capacity, signature='(n),(m)->(m)')

        refr_heat_capacities = np.vectorize(
            self.refrigerant_mix.mix_heat_capacity)
        
        sol = solve_bvp(
            odesystem, bc, z, in_guess, tol=tol, 
            max_nodes=max_nodes,verbose=verbose
        )
        
        return sol