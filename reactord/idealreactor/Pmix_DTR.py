from ..ReactorBase import ReactorBase
from scipy.integrate import solve_bvp
import numpy as np



class Pmix_DTR(ReactorBase):
    def __init__(self, mix, kinetic, t_operation, p_operation,
                time_minmax, f_in, f_out, t_in, 
                n_initial, t_initial, p_initial,
                n_final, t_final, p_final,
                q=0.0
    ):
        ReactorBase.__init__(self, mix, kinetic, t_operation, p_operation,
                            time_minmax=time_minmax,
                            f_in=f_in, f_out=f_out,
                            t_in=t_in, 
                            n_initial=n_initial, t_initial=t_initial, 
                            p_initial=p_initial, n_final=n_final, 
                            t_final=t_final, p_final=p_final,
                            q=0.0
        )

    def _grid_builder(self, N):
        time_array = np.linspace(self._time_minmax[0], 
                                self._time_minmax[1], 
                                N
        )
        return time_array

    def _border_condition_builder(self, ya, yb):
        bc = np.array([])
        for i,(n_in, n_fi) in enumerate(zip(self._n_initial, self._n_final)):
            if n_in == 'var':
                bc = np.append(bc, yb[i] - n_fi)
            else:
                bc = np.append(bc, ya[i] - n_in)

        if self._thermal_operation.lower() == 'isothermal':
            bc = np.append(bc, ya[-2] - self._t_initial)
        elif self._thermal_operation.lower() == 'non-isothermal':
            if self._t_initial == 'var':
                bc = np.append(bc, yb[-2] - self._t_final)
            else:
                bc = np.append(bc, ya[-2] - self._t_initial)

        if self._pressure_operation.lower() == 'isobaric':
            bc = np.append(bc, ya[-1] - self._p_initial)
        elif self._thermal_operation.lower() == 'non-isobaric':
            if self._t_initial == 'var':
                bc = np.append(bc, yb[-1] - self._p_final)
            else:
                bc = np.append(bc, ya[-1] - self._p_initial)
                
        return bc

    def _initial_guess_builder(self, N):
        n_comp = len(self._mix.substances)
        in_guess = np.zeros([n_comp+1,N])
        
        for i,(n_in, n_fi) in enumerate(zip(self._n_initial, self._n_final)):
            if n_in == 'var':
                in_guess[i,:] = np.full((N,1), n_fi)
            else:
                in_guess[i,:] = np.full((N,1), n_in)
        
        if self._thermal_operation.lower() == 'isothermal':
            in_guess[-1,:] = np.full((N,1), self._t_initial)
        elif self._thermal_operation.lower() == 'non-isothermal':
            if self._t_in == 'var':
                in_guess[-1,:] = np.full((N,1), self._t_final)
            else:
                in_guess[-1,:] = np.full((N,1), self._t_initial)
        return in_guess

    def _mass_balance(self,time,ns,*args):
        t,p,N = args
        n_comp = len(self._mix.substances)

        #Reaction volume calculation
        v = self._mix.volume(t,p)
        #kinetic evaluation on each of the N discretizations of z
        kinetics_eval =  np.array([self.kinetic.func(ns[:,i],t[i],p) 
                                  for i in range(N)
        ])
                
        #Mass balance of the n substances in mix
        dn_dtime = np.array([kinetic * v for kinetic in kinetics_eval])
        
        return np.vstack(dn_dtime)

    def _energy_balance(self,time,t,*args):
        ns,p,N = args
        
        v = self._mix.volume(t,p)

        if self._thermal_operation.lower() == 'isothermal':
            dT_dtime = 0.0
        elif self._thermal_operation.lower() == 'non-isothermal':
            DH = self._kinetic.DH
            
            ri = np.array([
                            self._kinetic.rs(ns[:,i],t[i],p[i]) 
                            for i in range(N)
            ])
            
            cp_mix = np.array([self._mix.mix_heat_capasity(ns[:,i],t[i],p[i]) 
                           for i in range(N)
            ])
            
            Q_time = self._q(time)
            
            dT_dtime = ((Q_time - np.dot(-DH, ri) * v) / (sum(ns) * cp_mix))
            
        return dT_dtime

    def _pressure_balance(self):
        pass

    def _particle_mass_balance(self):
        pass
    
    def _particle_energy_balance(self):
        pass
    
    def _result_report(self):
        pass

    def solve(self, N ,tol=0.001, max_nodes=1000, verbose=0):
        z_array = self._grid_builder(N)
        bc = self._border_condition_builder
        in_guess = self._initial_guess_builder(N)

        def odesystem(self,z,vars):
            ns = vars[0:len(vars)-1]
            t = vars[-1]
            #Mass balance
            dns_dtime = self._mass_balance(z,ns,[t, self._p, N])
            #Energy balance
            dt_dtime = self._energy_balance(z,t,ns,self._p)

            dvars_dz = np.append(dns_dtime, dt_dtime)
            return dvars_dz

        sol = solve_bvp(odesystem, bc, z_array, in_guess, 
                        tol=tol,
                        max_nodes=max_nodes,  
                        verbose=verbose)

        return z_array, sol.y