from ReactorBase import ReactorBase
from scipy.integrate import solve_bvp

class PFE(ReactorBase):
    def __init__(self, kinetic, f_ins, f_out, t_initial, t_final, p_initial,
                 p_final, kinetic_arg, operation, Q=None):
        
        ReactorBase.__init__(kinetic=kinetic, f_ins=f_ins, f_out=f_out,
        t_initial=t_initial, t_final=t_final, p_initial=p_initial, 
        p_final=p_final, kinetic_arg=kinetic_arg, operation=operation, Q=Q)  
                             