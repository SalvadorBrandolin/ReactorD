import numpy as np

class Mix:
    """ Set up of the initial conditions in the reactor.

    Parameters
    ----------
    
    component_list: List (strings) 
        The different compounds in the reactor
    phase: string
        Pashe in which the reaction take place. Can be "liquid" or "gas"
    moles: List [mol]
        Intial moles of each compound
    temp: float [K]
        Initial temperature
    volume: float [m^3]
        Reactor volume
    T_p: float [Pa]
        Total pressure in the reactor
    Par_p: List (float) [Pa]
        Partial pressures of all the compounds
    molar_frac: list [float] [adimensional]
        the molar fraction of each  in the reactor
    total_moles: float [moles]
        Total moles in the reactor
    concentrations: List (float) [moles]
        A list of the initial concentrations of the compounds in the reactor

    
    """
    

    def __init__(self, component_list, *, phase=None, moles=None, temp=None, volume=None, T_p=None,
                Par_p=None, molar_frac=None, total_moles=None, concentrations= None ):
        R= 8.31446261815324 # J/mol.K
        self.component_list= component_list
        self.phase= phase.lower()
        self.moles= moles        
        self.temp= temp
        self.volume= volume
        self.T_p= T_p
        self.Par_p= Par_p
        self.molar_frac= molar_frac
        self.total_moles= total_moles
        self.concentrations= concentrations
                      
        if self.concentrations == None:
            self.concentrations= np.zeros(len(self.component_list))
        if self.Par_p == None:
            self.Par_p=np.zeros(len(self.component_list))
        if self.moles == None:
            self.moles=np.zeros(len(self.component_list))
        
        #AGREGAR UNA COMPROBACION DE QUE LA SUMA DE Xi<=1

        #Concentrations given initial moles and volume:
        #if self.volume != None and self.moles != None:
        if self.volume !=None and ~any(self.moles):            
            for idx, component in enumerate (self.component_list):
                if self.volume <= 0:
                    print ("Volume must be greater than zero")
                    break
                self.concentrations[idx]= self.moles[idx] / self.volume
                #print(idx, component, self.concentrations[idx])

        #Concentrations given total moles, volume and molar fractions: 
        if self.total_moles != None and self.molar_frac != None and self.volume != None:
            ### SUMA DE FRACCIONES MOLARES sum(self.molar_frac)<=1.0
            for idx, component in enumerate (self.component_list):
                if self.volume <= 0:
                    print ("Volume must be different from zero")
                    break
                self.concentrations[idx]= self.total_moles * self.molar_frac[idx] / self.volume        
        
        #Partial pressures given total pressure and molar fractions:
        #We enter the IF when self.Par_p has any non-negative element
        if not(any(self.Par_p)) and self.molar_frac != None and self.T_p != None:                            
            for idx, component in enumerate (self.component_list):
                self.Par_p[idx]= self.molar_frac[idx] * self.T_p
            
        #Concentrations given partial pressures:
        #We enter the IF when self.Par_p has any non-negative element
        if (any(self.Par_p)) and self.temp!=None and self.T_p != None:            
            for idx, component in enumerate (self.component_list):
                self.concentrations[idx]= self.Par_p[idx] / (R*self.temp)