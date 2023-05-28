"""Module to calculate reaction enthalpy."""
from typing import TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .kinetic import Kinetic


def dh_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray[np.float64], float],
    pressure: Union[NDArray[np.float64], float],
) -> NDArray[np.float64]:
    """Return the user specified reaction enthalpies.

    Parameters
    ----------
    kinetic : Kinetic
        kinetic object
    temperature : Union[NDArray[np.float64], float]
        Temperature [K]
    pressure : Union[NDArray[np.float64], float]
        Pressure [Pa]

    Returns
    -------
    NDArray
        array of enthalpies of reactions set by user.
    """
    dh = np.full(
        (len(kinetic), np.size(temperature)),
        fill_value=kinetic._user_r_dhs[:, np.newaxis],
    )

    return dh


def dh_not_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray[np.float64], float],
    pressure: Union[NDArray[np.float64], float],
) -> NDArray[np.float64]:
    """Calculate the reaction enthalpies.

    Uses standard formation enthalpies and heat capacities to calculate the
    reaction enthalpies.

    Parameters
    ----------
    kinetic : Kinetic
        kinetic object
    temperature : Union[NDArray[np.float64], float]
        Temperature [K]
    pressure : Union[NDArray[np.float64], float]
        Pressure [Pa]

    Returns
    -------
    NDArray
        array of enthalpies of reactions set by methods in mix class.
    """
    corrs = kinetic.mix.formation_enthalpies_correction(temperature, pressure)

    if np.size(temperature) > 1:
        corrs = corrs.reshape((len(kinetic.mix), np.size(temperature)))

    dh = (
        kinetic.std_reaction_enthalpies
        + np.matmul(kinetic.stoichiometry, corrs).T
    ).T

    return dh
