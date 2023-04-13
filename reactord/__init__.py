"""ReactorD.

Chemical Reactor Simulations and Design.
"""
from . import flowreactors, idealreactor, mix, substance
from .kinetics import Kinetics
from .reactorbase import ReactorBase
from .substance import Substance
from .utils import vectorize

__all__ = [
    "flowreactors",
    "idealreactor",
    "mix",
    "substance",
    "Kinetics",
    "ReactorBase",
    "vectorize",
    "Substance",
]
