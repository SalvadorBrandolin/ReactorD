
class Kinetics:
    """Kinetics object class"

    Parameters
    ----------
    ri: list of all kinetics reactions laws in the system.
    stoichiometry: matrix of coefficients of reagents.
    Each row is a reaction that must start with the coefficient
    goverment react in the kinetics law. Columns are the compounds in the reaction.
    """


    def __init__(self, ri, stoichiometry):
        self.ri=ri
        self.stoichiometry= stoichiometry

    def 
