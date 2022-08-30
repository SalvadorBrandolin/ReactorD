import numpy as np
#import FinsMixer
from pde import PDEBase

class ReactorBase(PDEBase):
    def __init__(self, kinetic, f_ins=None, f_out=None, n_initial=None,
                n_final=None, t_initial=None, t_final=None, p_initial=None,
                p_final=None, kinetic_arg=None, operation=None,Q=None  
        ):

        self.n_substances = len(kinetic)
        self.substances = kinetic.substances
        self.kinetic = kinetic
        #self.f_in = flow_mixer(f_ins) #PROGRAMAR UN MEZCLADOR DE ENTRADAS
        self.f_out = f_out
        self.n_initial = n_initial
        self.n_final = n_final
        self.t_initial = t_initial
        self.t_final = t_final
        self.p_initial = p_initial
        self.p_final = p_final
        self.kinetic_arg = kinetic_arg
        self.operation = operation
        
        """
        self.in_flows = []
        if f_in != None:
            self.in_flows.append(flow for flow in f_in.molar_flows)
        else:
            raise Exception("Input inlet molar flow currents")

        self.out_flows = []
        if f_out != None:
            self.out_flows.append(flow for flow in self.f_out.molar_flows)
        else:
            raise Exception("Input output molar flow currents")

        """

    def grid_generator(self):
        return None

    def initial_state_generator(self):
        return None

    def mass_balance(self):
        return None

    def energy_balance(self):
        return None

    def pressure_balance(self):
        return None

    def solve(self):
        return None
