import numpy as np


def linear(
    mixture,
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
    """
