from pde import PDEBase, 

"""
Vamos a resolver la ODE dx/dt = k*(1-x)
con k = 10 y como condiciones inicial: x=0 a t = 0 
"""
class Reaction(pde.PDEBase):
    def __init__(self, k, bc):
        self.k = k

    def evolution_rate(self, state, t=0):
        x = state #Del estado tomo la variable x

        x_t = self.k * (1 - x)

        return pde.FieldCollection([x_t])


grid = pde.grids