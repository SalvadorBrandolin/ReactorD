# -*- coding: utf-8 -*-
from re import T
import numpy as np
import thermo.chemical

class Substance:
    """"Substance" object defining class

    Parameters
    ----------

    name : string
         Name of the substance.
    mw : float
         Molecular weight of the subtance [g / mol].
    """

    def __init__(self, id=None, name=None, mw=None, Tc=None, Pc=None, w=None,
                VolumeLiquid_func):

        if id == None:
          pass
        else:
            self.thermo_chemical_object = thermo.chemical.Chemical(ID = id)

            #Pure compound properties:
            self.name = self.thermo_chemical_object.name
            self.mw = self.thermo_chemical_object.MW
            self.Tc = self.thermo_chemical_object.Tc
            self.Pc = self.thermo_chemical_object.Pc
            self.omega = self.thermo_chemical_object.omega

            #Temperature dependent properties calculation functions:
            VolumeLiquid_func = self.thermo_chemical_object.VolumeLiquid

    def VolumeLiquid(self, VolumeLiquid_func, T):
        return VolumeLiquid_func(T)
