import numpy as np


class Mix:
    """Mixture object generator class.

        Parameters
        ----------
        substance_list : ndarray or list[Substance objects]
            list or array of Substance objects.            
        phase : string
            string that indicates the phase nature of the mixture. The
            avaliable options are: 'liquid', 'gas'.
        """
    
    def __init__(self, substance_list, phase):
        self.substances = substance_list
        self.phase = phase.lower()
        
        # Initialization of the heats of formation
        if self.phase == 'liquid':
            self.h_formations = [
                substance.h_formation for substance in self.substances
            ]
        elif self.phase == 'gas':
                self.h_formations = [
                substance.h_formation_ig for substance in self.substances
            ]
        else:
            raise ValueError(f'{self.phase} is not a supported phase') 
        
    def concentrations(self, moles, temperature, pressure):
        """Concentrations of the mixtures substances at the given moles 
        of each compound, temperature and pressure.
        
        Parameters:
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]
                
        Returns
        -------
        ndarray [float]
            ndarray that contains the concentrations of the mixture's 
            substances [mol/m^3]
        """

        zi = self.mol_fracations(moles)
        
        if self.phase == 'liquid':
            molar_volumes = np.array(
                [substance.volume_liquid(temperature, pressure) 
                for substance in self.substances]
            )                        
            
        elif self.phase == 'gas':
            molar_volumes = np.array(
                [substance.volume_gas(temperature, pressure) 
                for substance in self.substances]
            )
            
        total_molar_vol = np.dot(molar_volumes, zi)
        concentrations = np.divide(zi, total_molar_vol) #moles/m^3     
        return concentrations
 
    def volume(self, moles, temperature, pressure):
        """Method that returns the volume of the mixture.

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]

        Returns
        -------
        float
            volume of the mixture [m^3]
        """
        if self.phase == 'liquid':
            pure_volumes = np.array(
                [substance.volume_liquid(temperature, pressure) 
                for substance in self.substances]
            )
            return np.dot(pure_volumes, moles)

        if self.phase == 'gas':
            pure_volumes = np.array(
                [substance.volume_gas(temperature, pressure) 
                for substance in self.substances]
            )
            return np.dot(pure_volumes, moles)
        
    def mix_heat_capacity(self, moles, temperature, pressure): 
        """Method that returns the heat capacity of the mixture.

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]

        Returns
        -------
        float
            heat capacity of the mixture [j/mol/K)] 
        """

        zi = self.mol_fracations(moles)

        if self.phase == 'liquid':
            pure_cp = np.array(
                [substance.heat_capacity_liquid(temperature) 
                for substance in self.substances]
            )
            mix_cp = np.dot(zi, pure_cp)
            return mix_cp

        elif self.phase == 'gas':
            pure_cp = np.array([
                                substance.heat_capacity_gas(temperature) 
                                for substance in self.substances]
            )
            mix_cp = np.dot(zi, pure_cp)
            return mix_cp

    def mol_fracations(self, moles):
        """method that calculates the molar fractions of the mixture

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance

        Returns
        -------
        ndarray
            array that contains the molar fractions of mixture's 
            substances
        """
        total_moles = np.sum(moles, axis=0)
        zi = np.divide(moles, total_moles)
        return zi

    def partial_pressures(self, moles, temperature, pressure):
        """method that calculates the partial pressures of the mixture

        Parameters
        ----------
        moles: ndarray or list [float]
            moles of each substance
        temperature: float
            Temperature [K]
        pressure: float
           Total Pressure [Pa]

        Returns
        -------
        ndarray
            array that contains the partial pressures of mixture's 
            substances
        """
        zi = self.mol_fracations(moles)
        partial_pressures= np.multiply(zi, pressure)
        return partial_pressures

    def partial_P_2_conc (self, partial_pressures, T):
        R= 8.31446261815324 # J/mol.K
        self.partial_pressures= np.array(partial_pressures)
        conc= self.partial_pressures /(R*T) # mol/m^3
        return conc

    def __len__(self):
        return len(self.substances)
    
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