import numpy as np
from numpy.typing import NDArray

def ergun_pressure_balance(pfr) -> NDArray:
        mix_molecular_weight = pfr.mix.mix_molecular_weight(pfr.mole_fractions)
        total_molar_flow = np.sum(pfr.molar_flow, axis=0)
        total_mass_flow = total_molar_flow * mix_molecular_weight / 1000

        mass_density = pfr.mix.mass_density(
            pfr.mole_fractions, pfr.temperature, pfr.pressure
        )
        mix_viscosity = pfr.mix.mixture_viscosity(
            pfr.mole_fractions, pfr.temperature, pfr.pressure
        )

        f = total_mass_flow / pfr.transversal_area
        rho = mass_density
        phi = pfr.packed_bed_porosity
        dp = pfr.packed_bed_particle_diameter
        mu = mix_viscosity

        dp_dz = (
            -f
            / (rho * dp)
            * (1 - phi)
            / phi**3
            * (150 * (1 - phi) * mu / dp + 1.75 * f)
        )
        return dp_dz