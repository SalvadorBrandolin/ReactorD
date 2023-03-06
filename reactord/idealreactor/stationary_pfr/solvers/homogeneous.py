

def homogeneous_solver(
        self, grid_size=1000, tol=0.001, max_nodes=1000, verbose=0
    ) -> List[float]:
        """Simulate the homogeneous reactor.

        Parameters
        ----------
        grid_size : int, optional
            _description_, by default 1000
        tol : float, optional
            tolerance argument for scipy.solve_bvp , by default 0.001
        max_nodes : int, optional
            max_nodes argument for scipy.solve_bvp, by default 1000
        verbose : int, optional
            verbose argument for scipy.solve_bvp. Options 0 (False) or 1
            (True), by default 0

        Returns
        -------
        List[float]
            Solution matrix.
        """

        def odesystem(
            length_coordinate: List[float], variables: List[float]
        ) -> List[float]:
            """ODE system for scipy.solve_bvp.

            Parameters
            ----------
            length_coordinate : List[float]
                Reactor's length coordinate. [m]
            variables : List[float]
                Dependent variables on a two dimensional array in the
                order:
                    variables: [
                        [flux_1(z1), flux_1(z2), ..., flux_1(zn)],
                        [flux_2(z1), flux_2(z2), ..., flux_2(zn)],
                        .
                        .
                        .
                        [flux_m(z1), flux_m(z2), ..., flux_m(zn)],
                        [   T(z1)  ,   T(z2)   , ...,   T(zn)   ],
                        [   P(z1)  ,   P(z2)   , ...,   P(zn)   ],
                        [ Tref(z1) ,  Tref(z2) , ...,  Tref(zn) ],
                    ]
                where:
                    flux_i = molar flux of substance i
                    z_j = length_coordinate j
                    T = reactor's temperature
                    P = reactor's pressure
                    Tref = refrigerant's temperature

            Returns
            -------
            List[float]
                Derivatives matrix at each length_coordinate.
            """
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

        grid = self._grid_builder(grid_size)

        bc, initial_guess = self._border_cond_and_initial_guesses(grid_size)

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