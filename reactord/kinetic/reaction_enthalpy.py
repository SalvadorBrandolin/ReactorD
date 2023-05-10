from typing import TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .kinetic import Kinetic


def dh_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray, float],
    pressure: Union[NDArray, float],
) -> NDArray:
    dh = np.full(
        (len(kinetic), np.size(temperature)),
        fill_value=kinetic._user_r_dhs[:, np.newaxis],
    )

    return dh


def dh_not_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray, float],
    pressure: Union[NDArray, float],
) -> NDArray:
    corrs = kinetic.mix.formation_enthalpies_correction(temperature, pressure)

    if np.size(temperature) > 1:
        corrs = corrs.reshape((len(kinetic.mix), np.size(temperature)))

    dh = (
        kinetic.std_reaction_enthalpies
        + np.matmul(kinetic.stoichiometry, corrs).T
    ).T

    return dh
