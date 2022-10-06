import numpy as np
from abc import ABCMeta, abstractmethod
from Substance import Substance
from thermo.chemical import Chemical


class Abstract_Mix(metaclass = ABCMeta):

    """Mixture object abstract class.

        Parameters
        ----------
        substance_list : ndarray or list[Substance objects]
            list or array of Substance objects."""
         
    def __init__(self, list_of_substances):
        self.substances = list_of_substances
        self.enthalpies_formation = self.enthalpies_formation_builder() 
        
    @abstractmethod
    def enthalpies_formation_builder(self):
        pass

    @abstractmethod
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
        pass

    @abstractmethod
    def volume(self):
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
        pass

    @abstractmethod
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
        pass

# Other methods (Inhereted but not implemented in subclasses)
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
    
    def __len__(self):
        return len(self.substances)
    
    def __str__ (self):
        string=(f"The mixture contains the following" 
                f" {len(self.substances)} components:\n")
        for i,substance in enumerate (self.substances):
            string = string + substance.name.capitalize() + "\n"     
        return string


class Liquid_Mix(Abstract_Mix):

    def enthalpies_formation_builder(self):
        enthalpies_formation = [
                substance.h_formation for substance in self.substances
            ]
        return enthalpies_formation

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)      
        molar_volumes = np.array(
        [substance.volume_liquid(temperature, pressure) 
        for substance in self.substances]
        )                        

        total_molar_vol = np.dot(zi,molar_volumes)
        concentrations = np.divide(zi, total_molar_vol) #moles/m^3     
        return concentrations

    def volume(self, moles, temperature, pressure):  
        pure_volumes = np.array(
            [substance.volume_liquid(temperature, pressure) 
            for substance in self.substances]
        )
        return np.dot(pure_volumes, moles)        

    def mix_heat_capacity(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_cp = np.array(
            [substance.heat_capacity_liquid(temperature) 
            for substance in self.substances]
        )
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp


class IdealGas_Mix(Abstract_Mix):

    def enthalpies_formation_builder(self):
        enthalpies_formation = [
                substance.h_formation_ig for substance in self.substances
            ]
        return enthalpies_formation
    
    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)      
        
        molar_volumes = np.array(
        [substance.volume_gas(temperature, pressure) 
        for substance in self.substances]
        )

        total_molar_vol = np.dot(zi,molar_volumes)
        concentrations = np.divide(zi, total_molar_vol) #moles/m^3     
        return concentrations

    def volume(self, moles, temperature, pressure):
        pure_volumes = np.array(
            [substance.volume_gas(temperature, pressure) 
            for substance in self.substances]
        )
        return np.dot(pure_volumes, moles) 

    def mix_heat_capacity(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_cp = np.array([
                            substance.heat_capacity_gas(temperature) 
                            for substance in self.substances]
        )
        mix_cp = np.dot(zi, pure_cp)
        return mix_cp

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

    def partial_P_2_conc (self, partial_pressures, temperature):
        R= 8.31446261815324 # J/mol.K
        self.partial_pressures= np.array(partial_pressures)
        conc= self.partial_pressures /(R*temperature) # mol/m^3
        return conc