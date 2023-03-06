import numpy as np


def herning_zipperer(
    mixture,
    molar_fractions: np.ndarray,
    temperature: float,
    pressure: float,
) -> float:
    r"""Calculate the mix's viscosity with the Herning-Zipperer mixing rule.

    Multiple mixture compositions can be specified by a moles matrix. Each
    column of the matrix represents each mixture and each row represents each
    substance's mole fractions. Also with temperature and pressure vectors
    following de NumPy broadcasting rules.

    Ref: Herning, F.; Zipperer, L. (1936). "German: Beitrag zur Berechnung
    der Zähigkeit Technischer Gasgemische aus den Zähigkeitswerten der
    Einzelbestandteile; English: Calculation of the Viscosity of Technical Gas
    Mixtures from Viscosity of the individual Gases". Das Gas- und Wasserfach.
    79 (1936): 49–54 and 69–73.

    .. math::
           ln({\mu_{mix}}_{(T,P)}) =
           \frac {\sum z_i \mu_i_{(T,P)} \sqrt(M_i)} {\sum z_i \sqrt(M_i)}

    | :math:`\mu_{mix}`: mix's viscosity.
    | :math:`z_i`: mole fraction of the mix's :math:`i`-th substance.
    | :math:`\mu_i`: viscosity of the mix's :math:`i`-th substance.
    | :math:`M_i`: molecular weight of the mix's :math:`i`-th substance.
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
        molecular_weight defined on each mix's Substance.
    """

    # Take pure molecular weights
    molecular_weights = np.array(
        [substance.molecular_weight for substance in mixture.substances]
    )

    sqrt_molecular_weights = np.sqrt(molecular_weights)

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
    not_sum_denominator = np.multiply(molar_fractions, sqrt_molecular_weights)
    not_sum_numerator = np.multiply(not_sum_denominator, pure_viscosities)

    mixture_viscosities = np.divide(
        not_sum_numerator.sum(axis=0), not_sum_denominator.sum(axis=0)
    )

    return mixture_viscosities
