import numpy as np


class Mix:
    """ Set up of the initial conditions in the reactor."""
    
    def __init__(self, substance_list, phase):
        self.substances= substance_list #A list of objects
        self.phase= phase.lower()
        
    def concentrations(self, moles, T, P, setConc=False):
        """" Concentrations given moles of each compound, T and P.
        The molar volume of each compound is calculated and the
        total molar volume is calculated and used to return the
        concentrations.
        Volumes are supposed to be aditive.

        Parameters:
        moles: float
            moles of each substance
        T: float
            Temperature [K]
        P: float
           Total Pressure [Pa]
        setConc: boolean
            Allows to set concentrations directly. By default False
        """
        
        self.moles = np.array(moles)     
        zi= self.mol_frac(self.moles)
        Total_Molar_Vol = 0   # [m^3/mol]
        if self.phase == 'liquid':                        
            for i, substance in enumerate (self.substances):
                Total_Molar_Vol = Total_Molar_Vol + substance.volume_liquid(T, P) * zi[i]
        elif self.phase == 'gas':
            for i, substance in enumerate (self.substances):
                Total_Molar_Vol = Total_Molar_Vol + substance.volume_gas(T, P) * zi[i]

        elif setConc == True: # PARA SETEAR CONCENTRACIONES A PARTIR DE LA DENSIDAD
            pass              # Y ZI, O ALGUNA OTRA OPCION. QUEDA PENDIENTE.........
        
        conc= zi / Total_Molar_Vol    #concentrations in moles/m^3     
        return conc
 
    def volume(self, moles, T, P):
        if self.phase == 'liquid':
            pure_volumes = np.array([substance.volume_liquid(T,P) for substance
                                    in self.substances 
            ])
            return sum(pure_volumes)

        if self.phase == 'gas':
            pure_volumes = np.array([substance.volume_gas(T,P) for substance
                                    in self.substances 
            ])
            return sum(pure_volumes)
        
    def mix_heat_capacity(self, moles, T, set=None):
        zi = self.mol_frac(moles)
        if self.phase == 'liquid':
            pure_cp = np.array([
                                substance.heat_capacity_liquid(self, T) for 
                                substance in self.substances
            ])
            mix_cp = np.dot(zi, pure_cp)
            return mix_cp

        elif self.phase == 'gas':
            pure_cp = np.array([
                                substance.heat_capacity_gas(self, T) for 
                                substance in self.substances
            ])
            mix_cp = np.dot(zi, pure_cp)
            return mix_cp

    def mol_frac(self, moles):
        """Molar fractions calculator. 
        zi are the molar fractions"""
        self.moles= np.array(moles)
        self.total_moles= sum(self.moles)
        zi= self.moles / self.total_moles
        return zi

    def partial_P(self, moles, P):
        """Partial Pressure calculations using molar fractions and total pressure.
        zi are the molar fractions"""
        self.moles= np.array(moles)
        zi= self.mol_frac(self.moles)
        Pp= zi * P
        return Pp

    def partial_P_2_conc (self, Pp, T):
        R= 8.31446261815324 # J/mol.K
        self.Pp= np.array(Pp)
        conc= self.Pp /(R*T) # mol/m^3
        return conc
    
    def __str__ (self):
        string=(f"The mixture is in {self.phase} phase and " 
                f"contains the following {len(self.substances)} components:\n")
        for i,substance in enumerate (self.substances):
            string = string + substance.name.capitalize() + "\n"     
        return string
        """Observaciones: Quise hacer que devuelva un listado de concentraciones,
        fracciones molares y presiones parciales cuando hacemos el print del
        objeto, pero para eso tendría que calcular todo eso dentro de la funcion
        __str__. Como los argumentos de las funciones implicadas NO son atributos
        de la clase Mix, no lo puedo hacer (mejor dicho, no sé como hacerlo, si es
        que se puede) 
        """


      


        


        


                

          