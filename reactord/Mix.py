import numpy as np
from thermo.eos import R


class Mix:
    """ Set up of the initial conditions in the reactor."""
    
    def __init__(self, substance_list, phase):
        self.substances = substance_list #A list of objects
        self.phase = phase.lower()
        
    def concentrations(self, moles, temp, pressure, setConc = False):
        """ Concentrations given moles of each compound, T and P.
        The molar volume of each compound is calculated and the
        total molar volume is calculated and used to return the
        concentrations.
        Volumes are supposed to be aditive.

        Parameters:
        moles: float
            moles of each substance
        temp: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]
        setConc: boolean
            Allows to set concentrations directly. By default False
        """
        
        self.moles = np.array(moles)     
        zi = self.mol_frac(self.moles) # zi are the molar fractions
        total_molar_vol = 0   # [m^3/mol]

        if self.phase == 'liquid':
            molar_volumes = np.array(
                [substance.volume_liquid(temp, pressure) 
                for substance in self.substances]
            )                        
            
        elif self.phase == 'gas':
            molar_volumes = np.array(
                [substance.volume_gas(temp, pressure) 
                for substance in self.substances]
            )
            
        elif setConc == True: # PARA SETEAR CONCENTRACIONES A PARTIR DE LA DENSIDAD
            pass              # Y ZI, O ALGUNA OTRA OPCION. QUEDA PENDIENTE.........
        
        total_molar_vol = np.dot(molar_volumes, zi)
        conc = zi / total_molar_vol    #concentrations in moles/m^3     
        return conc
 
    def volume(self, moles, temp, pressure):
        if self.phase == 'liquid':
            pure_volumes = np.array([substance.volume_liquid(temp, pressure) for substance
                                    in self.substances 
            ])
            
        if self.phase == 'gas':
            pure_volumes = np.array([substance.volume_gas(temp, pressure) for substance
                                    in self.substances 
            ])
            return np.dot(pure_volumes, moles)
        
    def mix_heat_capacity(self, moles, temp, set=None):
        zi = self.mol_frac(moles)
        if self.phase == 'liquid':            
            pure_cp = np.array([
                                substance.heat_capacity_liquid(temp) for 
                                substance in self.substances
            ])            
        elif self.phase == 'gas':
            pure_cp = np.array([
                                substance.heat_capacity_gas(temp) for 
                                substance in self.substances
            ])
            
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp

    def mol_frac(self, moles):
        """Molar fractions calculator. 
        zi is a np array with the molar fractions"""
        self.moles = np.array(moles)
        self.total_moles = sum(self.moles)
        zi = self.moles / self.total_moles
        return zi

    def partial_p(self, moles, pressure):
        """Partial pressure calculations using molar fractions and total pressure.
        zi is a np array with the molar fractions"""
        self.moles = np.array(moles)
        zi = self.mol_frac(self.moles)
        par_p = zi * pressure
        return par_p

    def partial_p2conc (self, par_p, temp):
        self.par_p = np.array(par_p)
        conc = self.par_p / (R*temp) # mol/m^3
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


      


        


        


                

          