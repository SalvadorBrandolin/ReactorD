from numpy.typing import NDArray

from scipy.integrate import solve_bvp

def homogeneous_bvp_solver(
        pfr, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0
    ) -> NDArray:

        def odesystem(
            length_coordinate: NDArray, variables: NDArray
        ) -> NDArray:
            n_subs = len(self.mix)

            molar_fluxes = np.transpose(variables[0:n_subs, :])
            temperature = variables[-3, :]
            pressure = variables[-2, :]
            refrigerant_temperature = variables[-1, :]

            mass_derivatives = self._mass_balance_func(
                length_coordinate, molar_fluxes, temperature, pressure
            ).T

            temperature_derivatives = self._energy_balance_func(
                length_coordinate,
                molar_fluxes,
                temperature,
                pressure,
                refrigerant_temperature,
            ).T

            pressure_derivatives = self._pressure_balance_func(
                length_coordinate, molar_fluxes, temperature, pressure
            ).T

            refrigerant_temperature_derivatives = (
                self._refrigerant_energy_balance(
                    length_coordinate, temperature, refrigerant_temperature
                )
            ).T

            return np.vstack(
                (
                    mass_derivatives,
                    temperature_derivatives,
                    pressure_derivatives,
                    refrigerant_temperature_derivatives,
                )
            )

        grid = pfr._grid_builder(grid_size)

        bc, initial_guess = pfr._border_cond_and_initial_guesses(grid_size)

        simulation = solve_bvp(
            odesystem,
            bc,
            grid,
            initial_guess,
            tol=tol,
            max_nodes=max_nodes,
            verbose=verbose,
        )

        return simulation