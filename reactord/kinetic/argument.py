import numpy as np


class KineticArgument:
    def __init__(self, names) -> None:
        
        # Dictionary with substance indexing
        self.ref_dict = {}
        for idx, name in enumerate(names):
            self.ref_dict[name] = idx
            
        self.values = np.full(np.size(names), None)
        
    def __getitem__(self, key: str):
        try:
            return self.values[self.ref_dict[key]]
        except KeyError:
            raise KeyError(
                f"There is no substance named '{key}' in the Kinetic object."
                " Please, check if the was some typing error in the"
                " substances' names in the reaction rate functions."
            )