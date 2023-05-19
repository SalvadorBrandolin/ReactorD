from typing import List

import numpy as np


class CompositionalArgument:
    def __init__(self, names: List[str]) -> None:
        # Dictionary with substances' index
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
                " Please, check if it was some typing error in the"
                " substances' names in the reaction rate functions."
            )
