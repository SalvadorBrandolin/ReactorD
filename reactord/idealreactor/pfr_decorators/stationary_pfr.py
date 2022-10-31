from reactord.base_decorator import BaseDecorator
from reactord.idealreactor import PFR

class StationaryPFR(BaseDecorator):
    def __init__(self, reactor: PFR) -> None:
        pass

pfr = PFR(1,2,3,4,5,6)

estacionario = StationaryPFR(pfr)

estacionario.

