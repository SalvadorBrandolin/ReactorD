"""Linear mixing rule."""
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from reactord.mix.abstract_mix import AbstractMix

def linear(
    mixture: "AbstractMix",
    molar_fractions: np.ndarray,
    temperature: float,
    pressure: float,
) -> float:
    r"""Calculate the mixture's viscosity with the linear mixing rule.

    Multiple mixture compositions can be specified by a moles matrix. Each
    column of the matrix represents each mixture and each row represents each
    substance's mole fractions. Also with temperature and pressure vectors
    following de NumPy broadcasting rules.

    .. math::
           {\mu_{mix}}_{(T,P)} = \sum z_i {\mu_i}_{(T,P)}

    | :math:`\mu_{mix}`: mix's viscosity.
    | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
    | :math:`\mu_i`: viscosity of the mix's :math:`i`-th substance.
    | :math:`T`: temperature.
    | :math:`P`: pressure.

    Parameters
    ----------
    mixture : AbstractMix
        Mixture object.
    mole_fractions : np.ndarray [float]
        Mole fractions of each substance specified in the same order as the
        mix's substances order.
    temperature: float
        Temperature. [K]
    pressure: float
        Pressure. [Pa]

    Returns
    -------
    float
        Mix's viscosity. [Pa s]


    Requires
    --------
        Pure viscosity (liquid or gas) defined on each mix's Substance.
    """
    # Take pure viscosities
    if mixture.phase_nature == "liquid":
        pure_viscosities = np.array(
            [
                substance.viscosity_liquid(temperature, pressure)
                for substance in mixture.substances
            ]
        )

    elif mixture.phase_nature == "gas":
        pure_viscosities = np.array(
            [
                substance.viscosity_gas(temperature, pressure)
                for substance in mixture.substances
            ]
        )

    else:
        raise ValueError(f"{mixture.phase_nature} is not a valid phase nature")

    # Mixing rule
    mixture_viscosities = np.multiply(pure_viscosities, molar_fractions).sum(
        axis=0
    )

    return mixture_viscosities
