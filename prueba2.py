#%%

from pde import (PDEBase, FieldCollection, CartesianGrid, ScalarField, 
                 ImplicitSolver, Controller
                )
from pde.tools.numba import jit
from numba import njit
import matplotlib.pyplot as plt
from pde import config

config['numba.debug'] = True
"""
Vamos a resolver la ODE dCa/dt = - Fv * dCa/dV - k * Ca
con k = 10 y como condiciones inicial: x=0 a t = 0. Con 0 <= t <= 10  
"""
# =============================================================================
#  Clases
# -----------------------------------------------------------------------------

class ReactiveSystem:
    def __init__(self, ks, inflows, reactions):
        self.ks = ks
        self.inflows = inflows
        self.reactions = reactions
        self.net_flow = sum(inflows)

    @jit
    def react(self, state):
        return [
            r(self.ks, state) for r in self.reactions
        ]
    
    def __len__(self):        
        return len(self.reactions)


class ReactorSolver(PDEBase):
    def __init__(
        self, reactive_system, V, bcs,
        n_points=10, nt=10, t_end=10, t=0, max_it=1000, tol=1e-5
    ):
        # Number of reactives (hence, balances)
        self.n_reactives = len(reactive_system)

        # System dimension
        grid = CartesianGrid([[0, V]], n_points)
        states = [
            ScalarField(grid, data=0)
            for i in range(self.n_reactives)
        ]
        self.state = FieldCollection(states)

        # Border conditions
        self.bcs = bcs

        # Reactive system
        self.reactive_system = reactive_system

        # Time ranges
        self.nt = nt
        self.t_end = t_end
        self.dt = t_end/nt

        # Solver tolerance
        self.max_iter = max_it
        self.tol = tol       
    
    def evolution_rate(self, state, t=0):
        # Agregué el índice cero en el gradiente, hay que indicar que es
        #  la derivada en X más allá de que sea la única variable

        # Calculate reaction rates
        n_reactions = len(self.reactive_system)
        ks = self.reactive_system.ks
        reactions = [
            react(ks, state) for react in self.reactive_system.reactions
            ]

        # Net flow (of course it shouldn't be like this)
        net_flow = self.reactive_system.net_flow

        # Whole balances
        system = [
            - net_flow * state[i].gradient(self.bcs[i])   # Flow balance
            + reactions[i]                                # Reactions
            for i in range(self.n_reactives)
        ]
        return FieldCollection(system)

    def _make_pde_rhs_numba(self, state):
        #attributes locally available

        #Net flow (this is bad implemented - supposed to be vol flow)
        net_flow = self.reactive_system.net_flow

        #Number of reactives
        n_reactives = len(self.reactive_system)

        #kinetic constans
        ks = self.reactive_system.ks

        #Reaction function of each compound

        reactions = [
            react(ks, state).data for react in self.reactive_system.reactions
            ]
        
        #Operands definition
        gradient_c = [
            state.grid.make_operator("gradient", bc=self.bcs[i]) for i 
            in range(n_reactives)
            ]

        @jit
        def pde_rhs(state_data, t=0):
            """ compiled helper function evaluating right hand side """
            state_grad = []

            for i in range(n_reactives):
                gradient = gradient_c[i]
                state_grad.append(gradient(state_data[i]))

            #Reaction function of each compound
          
            """system = [
                - net_flow * state_grad[i]   # Flow balance
                + reactions[i]               # Reactions
                for i in range(n_reactives)
                ]"""

            system = []
            for i in range(n_reactives):
                system.append(- net_flow * state_grad[i])
            return system

        return FieldCollection(pde_rhs)

    def resolve(self, numba = False):
        if numba == False:
            solver = ImplicitSolver(self, self.max_iter, self.tol, 
                                    backend='numpy')
            controller = Controller(solver, t_range=self.t_end, 
                                    tracker=None)
            return controller.run(self.state)
        else:
            solver = ImplicitSolver(self, self.max_iter, self.tol, 
                                    backend='numba')
            controller = Controller(solver, t_range=self.t_end, 
                                    tracker=None)
            return controller.run(self.state)

# =============================================================================
# =============================================================================
#  Planteo de problema
# -----------------------------------------------------------------------------

# Variables generales
kr = 1
t0 = 0
tf = 10
V = 1
Fv = 1
Ca0 = 1

# Definicion de sistema reactivo
rs = ReactiveSystem(
    ks=[kr],
    inflows=[Fv, 0],
    reactions=[
        lambda ks, concentrations: -ks[0] * concentrations[0],
        lambda ks, concentrations:  ks[0] * concentrations[0],        
    ]
)


# Condiciones de borde
##  - x0 = CaO
##  - xf = dCa/dx = 0
bc_a = [
    {"value": Ca0},
    {"derivative": 0}
    ] 
bc_b = [
    {"value": 0},
    {"derivative": 0}
]

# Reactor
reactor = ReactorSolver(
    reactive_system=rs,
    V=V,
    bcs=[bc_a, bc_b],
    n_points = 100,
    t_end=10,
    nt = 1000,
)
# ==============================================================================

results = reactor.resolve(numba=True)


plt.plot(results[0].data)
plt.plot(results[1].data)
plt.show()




# %%
