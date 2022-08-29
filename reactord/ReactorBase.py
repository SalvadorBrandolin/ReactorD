import numpy as np
#import FinsMixer
from pde import PDEBase

class PDEReactorBase(PDEBase):
    def __init__(self, kinetic, x_spans=None, Ns=None, N0=None, Nf=None, 
        f_ins=None, f_out=None, kinetic_arg=None, operation=None, T=None, 
        P=None,  
        ):

        self.n_substances = len(kinetic)
        self.substances = kinetic.substances
        self.kinetic = kinetic
        self.operation = operation
        self.T = T
        self.P = P
        self.N = Ns
        self.x_span = x_spans
        #f_in = flow_mixer(f_ins) #PROGRAMAR UN MEZCLADOR DE ENTRADAS
        f_out = f_out
        
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
