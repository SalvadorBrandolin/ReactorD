class Mix:
    """
    """ 

    def __init__(self, component_list, moles, phase):
        self.phase= phase.lower()
        self.component_list= component_list
        self.moles= moles

        if self.phase== 'liquid':
            for idx, components in enumerate (self.component_list):
                self.concentrations[idx]=
