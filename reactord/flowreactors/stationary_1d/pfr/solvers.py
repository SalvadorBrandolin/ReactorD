from typing import TYPE_CHECKING

import numpy as np

import pandas as pd

from scipy.integrate import solve_ivp, solve_bvp


if TYPE_CHECKING:
    from reactord.flowreactors.stationary_1d.pfr.pfr import PFR


def simulate_ivp_problem(reactor: "PFR", options:dict = {}):
    reactor.grid_size = reactor._initial_grid_size
    reactor.initial_profile_builder()
    reactor.border_conditions_builder()
    
    if options.get("t_eval") is None:
        t_eval = options.get("t_eval")
    else:
        t_eval = np.linspace(0, reactor.reactor_length, reactor.grid_size)
        
    if options.get("method"):
        method = options.get("method")
    else:
        method = "RK45"
    
    reactor.ode_solution = solve_ivp(
    fun = reactor.evaluate_balances,
    t_span = np.array([0, reactor.reactor_length]),
    y0 = reactor.inlet_conditions,
    method = "RK45",
    t_eval = t_eval,
    vectorized=True
    )
    
    # Update profiles with solution
    reactor.z = reactor.ode_solution.t
    reactor.mass_balance.update_profile(reactor, reactor.ode_solution.y)
    reactor.energy_balance.update_profile(reactor, reactor.ode_solution.y)
    reactor.pressure_balance.update_profile(reactor, reactor.ode_solution.y)
    reactor.mole_fraction_profile = reactor.mix.mole_fractions(
        reactor.mass_profile
    )
    reactor.r_rates_profile = reactor.kinetic.evaluate(
        reactor.mole_fraction_profile,
        reactor.temperature_profile,
        reactor.pressure_profile,
    )

    # Build data frame
    result = np.vstack((reactor.z, reactor.ode_solution.y))
    z = np.array(["z"])
    names = reactor.mix.names
    if reactor.refrigerant_temperature_profile is not None:
        last = np.array(
            ["temperature", "refrigerant_temperature", "pressure"]
        )
    else:
        last = np.array(["temperature", "pressure"])

    columns = np.concatenate((z, names, last))
    reactor.sim_df = pd.DataFrame(result.T, columns=columns, index=None)


def simulate_bvp_problem(reactor: "PFR", options:dict = {}):
        reactor.grid_size = reactor._initial_grid_size
        reactor.initial_profile_builder()
        reactor.border_conditions_builder()

        if options.get("tol"):
            tol = options.get("tol")
        else:
            tol = 1e-3
            
        if options.get("max_nodes"):
            max_nodes = options.get("max_nodes")
        else:
            max_nodes = 1000
            
        if options.get("verbose"):
            verbose = options.get("verbose")
        else:
            verbose = 0
            
        if options.get("bc_tol"):
            bc_tol = options.get("bc_tol")
        else:
            bc_tol = None

        # Simualate
        reactor.ode_solution = solve_bvp(
            fun=reactor.evaluate_balances,
            bc=reactor.border_conditions,
            x=reactor.z,
            y=reactor.initial_variables_profile,
            tol=tol,
            max_nodes=max_nodes,
            verbose=verbose,
            bc_tol=bc_tol,
        )

        # Update profiles with solution
        reactor.z = reactor.ode_solution.x
        reactor.mass_balance.update_profile(reactor, reactor.ode_solution.y)
        reactor.energy_balance.update_profile(reactor, reactor.ode_solution.y)
        reactor.pressure_balance.update_profile(reactor, reactor.ode_solution.y)
        reactor.mole_fraction_profile = reactor.mix.mole_fractions(
            reactor.mass_profile
        )
        reactor.r_rates_profile = reactor.kinetic.evaluate(
            reactor.mole_fraction_profile,
            reactor.temperature_profile,
            reactor.pressure_profile,
        )

        # Build data frame
        result = np.vstack((reactor.z, reactor.ode_solution.y))
        z = np.array(["z"])
        names = reactor.mix.names
        if reactor.refrigerant_temperature_profile is not None:
            last = np.array(
                ["temperature", "refrigerant_temperature", "pressure"]
            )
        else:
            last = np.array(["temperature", "pressure"])

        columns = np.concatenate((z, names, last))
        reactor.sim_df = pd.DataFrame(result.T, columns=columns, index=None)