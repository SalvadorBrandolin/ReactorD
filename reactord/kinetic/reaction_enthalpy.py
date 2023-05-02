from typing import TYPE_CHECKING, Union

from numpy.typing import NDArray

if TYPE_CHECKING:
    from .kinetic import Kinetic


def dh_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray, float],
    pressure: Union[NDArray, float],
) -> NDArray:
    return kinetic.r_dh


def dh_not_specified(
    kinetic: "Kinetic",
    temperature: Union[NDArray, float],
    pressure: Union[NDArray, float],
) -> NDArray:
    corrs = kinetic.mix.formation_enthalpies_correction(temperature, pressure)
    return kinetic.std_reaction_enthalpies + corrs
