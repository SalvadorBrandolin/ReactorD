"""ReactorD.

Chemical Reactor Simulations and Design.
"""
from . import flowreactors, idealreactor, kinetic, mix, substance
from .kinetic import Kinetic
from .reactorbase import ReactorBase
from .substance import Substance
from .utils import vectorize

__all__ = [
    "flowreactors",
    "idealreactor",
    "kinetic",
    "mix",
    "substance",
    "Kinetic",
    "ReactorBase",
    "vectorize",
    "Substance",
]
