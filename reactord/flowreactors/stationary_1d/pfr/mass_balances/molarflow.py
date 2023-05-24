"""PFR Molar flow module."""
from IPython.display import display

import numpy as np

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class MolarFlow:
    """PFR Molar Flow class.

    Parameters
    ----------
    molar_flows_in : dict, optional
        Molar flow in the inlet reactor section
    molar_flows_out : dict, optional
        Molar flow in the out reactor section
    """

    def __init__(
        self, molar_flows_in: dict = {}, molar_flows_out: dict = {}
    ) -> None:

        self.molar_flows_in = molar_flows_in
        self.molar_flows_out = molar_flows_out

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()))

    def initial_profile(self, reactor: PFR):
        """Set initial profile in PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            initial mass profile in all grid

        Raises
        ------
        ValueError
            Raise to guarantee one condition (inlet or outlet)
            for each substance
        ValueError
             Raise to guarantee border condition
        """
        n_substances = len(reactor.mix)
        initial_mass_profile = np.zeros((n_substances, reactor.grid_size))

        inlet_substances = list(self.molar_flows_in.keys())
        outlet_substances = list(self.molar_flows_out.keys())

        for idx, name in enumerate(reactor.mix.names):
            if name in inlet_substances and name in outlet_substances:
                raise ValueError(
                    "Mass balance error: both inlet and outlet border "
                    "conditions for the same substance not allowed. Check "
                    f"border conditions of the substance: {name}."
                )
            elif name in inlet_substances:
                initial_mass_profile[idx, :] = np.full(
                    reactor.grid_size, self.molar_flows_in[name]
                )
            elif name in outlet_substances:
                initial_mass_profile[idx, :] = np.full(
                    reactor.grid_size, self.molar_flows_out[name]
                )
            else:
                raise ValueError(
                    "Mass balance error: No border condition found for: "
                    f"{name}."
                )

        return initial_mass_profile

    # TODO
    def update_profile(self, reactor: PFR, variables):
        """Update profile."""
        reactor.mass_profile = variables[0 : reactor.subs_n, :]  # noqa

    def border_conditions(self, reactor: PFR):
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            array with molar flow in the inlet and outlet reactor
        """
        flows_in = np.array([])
        flows_out = np.array([])

        for substance_name in reactor.mix.names:
            flows_in = np.append(
                flows_in, self.molar_flows_in.get(substance_name)
            )
            flows_out = np.append(
                flows_out, self.molar_flows_out.get(substance_name)
            )

        return flows_in, flows_out

    def evaluate_balance(self, reactor: PFR):
        """Evaluate mass balance.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        ndarray
            Reaction rate of each substance on each reactor's z
        """
        # Reaction rate of each substance on each reactor's z
        ri = np.matmul(
            reactor.kinetic.stoichiometry.T, reactor.r_rates_profile
        )

        return ri * reactor.transversal_area

    def __repr__(self) -> str:
        """Represent equation of PFR mass balance."""
        return r"\frac{1}{a_t}\frac{dF_i}{dz}=r_i"
