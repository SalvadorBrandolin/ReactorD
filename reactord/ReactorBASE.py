from pde import PDEBase

class ReactorBASE(PDEBase):
    def __init__(self):
        #COSAS
        return
    
    def grid_state_gen(self):
        return

    def evolution_rate(self, state: FieldBase, t: float = 0) -> FieldBase:
        return super().evolution_rate(state, t)

    def solve(self):
        return


    