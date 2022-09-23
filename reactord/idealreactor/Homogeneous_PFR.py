from .. import ReactorBase
from scipy.integrate import solve_bvp
import numpy as np


class Homogeneous_PFR(ReactorBase):
    
    def __init__(
        self, mix, kinetic, t_operation,
        r_dims_minmax, 
        f_in,
        f_out,
        t_in, t_out,
        q,
        p,
        transversal_a=1,
        ):

        ReactorBase.__init__(
                            self, 
                            mix=mix, 
                            kinetic=kinetic, 
                            t_operation=t_operation, 
                            p_operation='isobaric',
                            r_dims_minmax=r_dims_minmax, 
                            f_in=f_in,
                            f_out=f_out,
                            t_in=t_in, t_out=t_out,
                            q=q)
        
        self._transversal_a = transversal_a
        self._p = p

    def _grid_builder(self, N):
        dim_array = np.linspace(self._r_dims_minmax[0], 
                                self._r_dims_minmax[1], 
                                N)
        return dim_array

    def _border_condition_builder(self, ya, yb):
        bc = np.array([])
        for i,(fin, fout) in enumerate(zip(self._f_in, self._f_out)):
            if fin == 'var':
                bc = np.append(bc, yb[i] - fout)
            else:
                bc = np.append(bc, ya[i] - fin)
        
        if self._thermal_operation.lower() == 'isothermal':
            bc = np.append(bc, ya[-1] - self._t_in)
        elif self._thermal_operation.lower() == 'non-isothermal':
            if self._t_in == 'var':
                bc = np.append(bc, yb[-1] - self._t_out)
            else:
                bc = np.append(bc, ya[-1] - self._t_in)
        return bc

    def _initial_guess_builder(self, N):
        n_comp = len(self._mix.substances)
        in_guess = np.zeros([n_comp+1,N])
        
        for i,(fin, fout) in enumerate(zip(self._f_in, self._f_out)):
            if fin is 'var':
                in_guess[i,:] = np.full((N,1), fout)
            else:
                in_guess[i,:] = np.full((N,1), fin)
        
        if self._thermal_operation.lower() == 'isothermal':
            in_guess[-1,:] = np.full((N,1), self._t_in)
        elif self._thermal_operation.lower() == 'non-isothermal':
            if self._t_in == 'var':
                in_guess[-1,:] = np.full((N,1), self._t_out)
            else:
                in_guess[-1,:] = np.full((N,1), self._t_in)
        return in_guess
                
    def _mass_balance(self,z,fs,*args):
        t,p,N = args

        n_comp = len(self._mix.substances)
        
        #kinetic evaluation on each of the N discretizations of z
        kinetics_eval =  np.array([self.kinetic.func(fs[:,i],t[i],p) 
                                  for i in range(N)])
        #Tranversal area:
        A = self._transversal_a
        
        #Mass balance of the n substances in mix
        dF_dz = np.array([A * kinetics_eval[i] for i in range(n_comp)])
        
        return np.vstack(dF_dz)

    def _energy_balance(self,z,t,*args):
        fs, N = args
        p = self._p
        
        A = self._transversal_a

        if self._thermal_operation.lower() == 'isothermal':
            dT_dz = 0.0
        elif self._thermal_operation.lower() == 'non-isothermal':
            DH = self._kinetic.DH
            
            ri = np.array([self._kinetic.rs(fs[:,i],t[i],p) for i in range(N)])
            
            cpm = np.array([self._mix.mix_heat_capasity(fs[:,i],t[i],p) 
                           for i in range(N)
                           ])
            
            Q_z = self._q(z)
            
            dT_dz = ((Q_z - np.dot(-DH, ri)) / (sum(fs) * cpm) * A)
            
        return dT_dz

    def _pressure_balance(self):
        pass

    def _particle_mass_balance(self):
        pass
    
    def _particle_energy_balance(self):
        pass

    def solve(self, N ,tol=0.001, max_nodes=1000, verbose=0):
        z_array = self._grid_builder(N)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(N)

        def odesystem(self,z,vars):
            fs = vars[0:len(vars)-1]
            t = vars[-1]
            #Mass balance
            dF_dz = self._mass_balance(z,fs,[t, self._p, N])
            #Energy balance
            dT_dz = self._energy_balance(z,t,fs,self._p)

            dvars_dz = np.append(dF_dz, dT_dz)
            return dvars_dz

        sol = solve_bvp(odesystem, bc, z_array, in_guess, 
                        tol=tol,
                        max_nodes=max_nodes,  
                        verbose=verbose)

        return z_array, sol.y