class ReadSystem:
    def __init__(
        self, reactive_system: dict, rate_argument: str = "concentration"
    ) -> None:
        self.r_names = list(reactive_system.keys())
        self.r_dicts = [reactive_system[name] for name in self.r_names]
        self.r_eqs = [rdict["eq"] for rdict in self.r_dicts]
        self.r_rates = [rdict["rates"] for rdict in self.r_dicts]
        self.r_dh = [rdict["DH"] for rdict in self.r_dicts]

    def __len__(self):
        return len(self.reaction_names)


d = {"esterificacion": {"eq": asd, "rate": asd, "DH": asd}}
