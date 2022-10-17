import numpy as np
from abc import ABCMeta, abstractmethod


class AbstractMix(metaclass=ABCMeta):
    """Mixture object abstract class.

    Parameters
    ----------
    substance_list : ndarray or list[Substance objects]
        list or array of Substance objects."""

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

    def __str__(self):
        string = (
            f"The mixture contains the following"
            f" {len(self.substances)} components:\n"
        )
        for i, substance in enumerate(self.substances):
            string = string + substance.name.capitalize() + "\n"
        return string
