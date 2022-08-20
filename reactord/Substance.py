# -*- coding: utf-8 -*-
import numpy as np

class Substance:
    """"Substance" object defining class

    Parameters
    ----------

    Name : string
         Name of the substance.
    Mw : float
         Molecular weight of the subtance [g / mol].
    """
    def __init__(self, Name = None, Mw = None, ):

        self.Name = Name
        self.Mw = Mw        