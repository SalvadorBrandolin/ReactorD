"""PFR module."""
import numpy as np
from numpy.typing import NDArray

import pandas as pd

from reactord.kinetic.kinetic import Kinetic

from scipy.integrate import solve_bvp


class PFR:
    """Piston flow reactor.

    No radial and no axial dispersion tubular reactor. Reaction on the tube
    side.

    Parameters
    ----------
    kinetic : Kinetic
        kinetic object.
    reactor_length : float
        Reactor length. [m]
    transversal_area : float
        Reactor transversal area. [m²]
    grid_size : int
        Number of discretization points of the reactor's length.
    mass_balance :
        Mass balance setting.
    energy_balance :
        Energy balance setting.
    pressure_balance :
        Pressure balance setting.
    """

    def __init__(
        self,
        kinetic: Kinetic,
        reactor_length: float,
        transversal_area: float,
        grid_size: int,
        mass_balance,
        energy_balance,
        pressure_balance,
    ) -> None:
        # =====================================================================
        # Core PFR information
        # =====================================================================
        self.kinetic = kinetic
        self.reactor_length = reactor_length
        self.transversal_area = transversal_area
        self._initial_grid_size = grid_size
        self.grid_size = grid_size
        self.subs_n = len(self.mix)
        self.reac_n = len(self.kinetic)

        self.tube_radius = np.sqrt(self.transversal_area / np.pi)

        # =====================================================================
        # Balances
        # =====================================================================
        self.mass_balance = mass_balance
        self.energy_balance = energy_balance
        self.pressure_balance = pressure_balance

        # =====================================================================
        # Variables grids and profiles
        # =====================================================================
        self.z = np.array([])
        self.mass_profile = np.array([])
        self.mole_fraction_profile = np.array([])
        self.temperature_profile = np.array([])
        self.refrigerant_temperature_profile = np.array([])
        self.pressure_profile = np.array([])
        self.r_rates_profile = np.array([])

    @property
    def mix(self):
        """Set mix of reaction."""
        return self.kinetic.mix

    def initial_profile_builder(self) -> None:
        """Construct initial profile of the reactor's dependent variables.

        Constructs the initial profile calling the initial_profile() method
        of each balance. This method sets the molar flow profile, reactor's
        temperature, refrigerant's temperature and reactor's pressure over the
        reactor's length.
        """
        self.z = np.linspace(0, self.reactor_length, self.grid_size)
        mass_profile = self.mass_balance.initial_profile(self)
        temperatures_profile = self.energy_balance.initial_profile(self)
        pressure_profile = self.pressure_balance.initial_profile(self)

        self.initial_variables_profile = np.vstack(
            (
                mass_profile,
                temperatures_profile,
                pressure_profile,
            )
        )

    def border_conditions_builder(self) -> None:
        """Extract the border conditions of each balance.

        Constructs the border conditions calling the border_conditions()
        method of each balance.
        """
        self.mass_bc = self.mass_balance.border_conditions(self)
        self.temperature_bc = self.energy_balance.border_conditions(self)
        self.pressure_bc = self.pressure_balance.border_conditions(self)

        # Inlet conditions:
        self.inlet_conditions = np.append(
            self.mass_bc[0], self.temperature_bc[0]
        )
        self.inlet_conditions = np.append(
            self.inlet_conditions, self.pressure_bc[0]
        )

        # Outlet conditions:
        self.outlet_conditions = np.append(
            self.mass_bc[1], self.temperature_bc[1]
        )
        self.outlet_conditions = np.append(
            self.outlet_conditions, self.pressure_bc[1]
        )

        # Index where border conditions are specified:
        self._in_index = np.argwhere(
            np.not_equal(self.inlet_conditions, None)
        ).ravel()
        self._out_index = np.argwhere(
            np.not_equal(self.outlet_conditions, None)
        ).ravel()

    def border_conditions(self, ya, yb) -> NDArray[np.float64]:
        """Border condition builder for scipy.solve_bvp.

        Parameters
        ----------
        ya : List[float]
            Border conditions on reactor's inlet.
        yb : List[float]
            Border conditions on reactor's outlet.

        Returns
        -------
        NDArray[np.float64]
            Border conditions for scipy.solve_bvp.
        """
        bc = np.array([])

        for idx in self._in_index:
            bc = np.append(bc, ya[idx] - self.inlet_conditions[idx])

        for idx in self._out_index:
            bc = np.append(bc, yb[idx] - self.outlet_conditions[idx])

        return bc

    def evaluate_balances(
        self, z: NDArray[np.float64], variables: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Evaluate mass, energy and pressure balance.

        Parameters
        ----------
        z: NDArray[np.float64]
            Reactor length
        Variables: NDArray[np.float64]
            Matrix of the molar flows, temperatures and pressure over the
            reactor's length.

        Returns
        -------
        NDArray[np.float64]
            Matrix of the gradients of molar flows, temperatures and pressure
            over the reactor's length.
        """
        self.z = z
        self.grid_size = np.size(z)
        self.mass_balance.update_profile(self, variables)
        self.energy_balance.update_profile(self, variables)
        self.pressure_balance.update_profile(self, variables)

        self.mole_fraction_profile = self.mix.mole_fractions(self.mass_profile)

        self.r_rates_profile = self.kinetic.evaluate(
            self.mole_fraction_profile,
            self.temperature_profile,
            self.pressure_profile,
        )

        mass_gradient = self.mass_balance.evaluate_balance(self)
        temperature_gradient = self.energy_balance.evaluate_balance(self)
        pressure_gradient = self.pressure_balance.evaluate_balance(self)

        gradients = np.vstack(
            (mass_gradient, temperature_gradient, pressure_gradient)
        )

        return gradients

    def simulate(
        self, tol=1e-3, max_nodes=1000, verbose=0, bc_tol=None
    ) -> None:
        """Simulate reactor.

        Parameters
        ----------
        tol : _type_, optional
            Scipy solve_bvp tol value, by default 1e-3
        max_nodes : int, optional
            Scipy solve_bvp max_nodes value, by default 1000
        verbose : int, optional
            Scipy solve_bvp verbose parameter, by default 0
        bc_tol : optional
            Scipy solve_bvp bc_tol value, by default None
        """
        self.grid_size = self._initial_grid_size
        self.initial_profile_builder()
        self.border_conditions_builder()

        # Simualate
        self.ode_solution = solve_bvp(
            fun=self.evaluate_balances,
            bc=self.border_conditions,
            x=self.z,
            y=self.initial_variables_profile,
            tol=tol,
            max_nodes=max_nodes,
            verbose=verbose,
            bc_tol=bc_tol,
        )

        # Update profiles with solution
        self.z = self.ode_solution.x
        self.mass_balance.update_profile(self, self.ode_solution.y)
        self.energy_balance.update_profile(self, self.ode_solution.y)
        self.pressure_balance.update_profile(self, self.ode_solution.y)
        self.mole_fraction_profile = self.mix.mole_fractions(self.mass_profile)
        self.r_rates_profile = self.kinetic.evaluate(
            self.mole_fraction_profile,
            self.temperature_profile,
            self.pressure_profile,
        )

        # Build data frame
        result = np.vstack((self.z, self.ode_solution.y))
        z = np.array(["z"])
        names = self.mix.names
        if self.refrigerant_temperature_profile is not None:
            last = np.array(
                ["temperature", "refrigerant_temperature", "pressure"]
            )
        else:
            last = np.array(["temperature", "pressure"])

        columns = np.concatenate((z, names, last))
        self.sim_df = pd.DataFrame(result.T, columns=columns, index=None)

    @property
    def irepr(self) -> None:
        """Represent for ipython."""
        self.kinetic.irepr
        print("Mass balance:")
        self.mass_balance.irepr
        print("Reactor and refrigerant energy balances:")
        self.energy_balance.irepr
        print("Pressure balance:")
        self.pressure_balance.irepr

    def __repr__(self) -> str:
        """Represent equation of PFR mass balance."""
        latex = (
            f"{self.kinetic.__repr__()}\n"
            f"{self.mass_balance.__repr__()}\n"
            f"{self.energy_balance.__repr__()[0]}\n"
            f"{self.energy_balance.__repr__()[1]}\n"
            f"{self.pressure_balance.__repr__()}\n"
        )
        return latex
