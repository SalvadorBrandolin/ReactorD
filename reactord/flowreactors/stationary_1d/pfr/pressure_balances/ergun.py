"""PFR non-isobaric pressure balance module."""
from IPython.display import display

import numpy as np
from numpy.typing import NDArray

from reactord.flowreactors.stationary_1d.pfr.pfr import PFR

from sympy import symbols


class Ergun:
    r"""PFR Ergun energy balance class.

    .. math::
        \frac{dP}{dz}=-\frac{G}{{\rho}D_p}\left(\frac{1-\phi}{\phi^3}
        \right)\left[\frac{150(1-\phi)\mu}{D_p}+1.75G\right]

    Parameters
    ----------
    pressure : dict
        pressure of each substance in the inlet or outlet position.
    porosity : float
        Packed bed porosity float between 0 and 1.
    particle_diameter : float
        Diameter of particle of packed bed. [m]

    Raises
    ------
    ValueError
        Only inlet or outlet border condition for pressure are allowed.
    """

    def __init__(
        self, pressure: dict, porosity: float, particle_diameter: float
    ) -> None:
        self._inlet_pressure = pressure.get("in")
        self._outlet_pressure = pressure.get("out")
        self.porosity = porosity
        self.particle_diameter = particle_diameter

        if self._inlet_pressure and self._outlet_pressure:
            raise ValueError(
                "Pressure balance error: Only inlet or outlet border condition"
                " specification allowed for the pressure."
            )

    @property
    def irepr(self):
        """Represent for ipython."""
        display(symbols(self.__repr__()))

    def initial_profile(self, reactor: PFR) -> NDArray[np.float64]:
        """Set initial pressure profile in non-isobaric PFR.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        NDArray[np.float64]
            initial pressure profile in all grid
        """
        if self._inlet_pressure:
            return np.full(reactor.grid_size, self._inlet_pressure)
        else:
            return np.full(reactor.grid_size, self._outlet_pressure)

    def update_profile(
        self, reactor: PFR, variables: NDArray[np.float64]
    ) -> None:
        """Update th reactor's presure profile.

        Parameters
        ----------
        reactor : PFR
            PFR object.
        variables : NDArray[np.float64]
            Variables of solve_bvp ode solver.
        """
        reactor.pressure_profile = variables[-1, :]

    def border_conditions(self, reactor: PFR) -> NDArray[np.float64]:
        """Set border conditions.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        NDArray[np.float64]
            array with pressure in the inlet and outlet reactor.
        """
        if self._inlet_pressure:
            return self._inlet_pressure, None
        else:
            return None, self._outlet_pressure

    def evaluate_balance(self, reactor: PFR) -> NDArray[np.float64]:
        """Evaluate pressure balance.

        Parameters
        ----------
        reactor : PFR
            PFR object

        Returns
        -------
        NDArray[np.float64]
            Pressure rate of each substance on each reactor's z.
        """
        phi = self.porosity
        dp = self.particle_diameter

        m_rho = reactor.mix.mass_density(
            reactor.mole_fraction_profile,
            reactor.temperature_profile,
            reactor.pressure_profile,
        )
        mu = reactor.mix.mix_viscosity(
            reactor.mole_fraction_profile,
            reactor.temperature_profile,
            reactor.pressure_profile,
        )
        g = (
            np.sum(reactor.mass_profile[:, 0])
            * reactor.mix.mix_molecular_weight(
                reactor.mole_fraction_profile[:, 0]
            )
            / 1000
            / reactor.transversal_area
        )

        pressure_gradient = (
            -g
            / m_rho
            / dp
            * (1 - phi)
            / phi**3
            * (150 * (1 - phi) * mu / dp + 1.75 * g)
        )
        return pressure_gradient

    def __repr__(self):
        """Represent equation of PFR non-isobaric pressure balance."""
        latex = (
            r"\frac{dP}{dz}=-\frac{G}{{\rho}D_p}\left(\frac{1-\phi}{\phi^3}"
            r"\right)\left[\frac{150(1-\phi)\mu}{D_p}+1.75G\right]"
        )
        return latex
