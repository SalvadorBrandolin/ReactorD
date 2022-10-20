from abc import ABCMeta, abstractmethod

import numpy as np


class AbstractMix(metaclass=ABCMeta):
    """Mixture object abstract class.

    Parameters
    ----------
    substance_list : ndarray or list[Substance objects]
        list or array of Substance objects."""

    @abstractmethod
    def concentrations(
        self, moles: list[float], temperature: float, pressure: float
    ):
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
        raise NotImplementedError()

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
        raise NotImplementedError()

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
        """
        raise NotImplementedError()

    @abstractmethod
    def formation_enthalpies(self):
        """Method that obtains/calculates the formation enthalpy of the
        mixture's substances.
        """
        raise NotImplementedError()

    @abstractmethod
    def formation_enthalpies_correction(
        self,
        temperature: float,
        pressure: float,
    ):
        """Method that correct the formation enthalpy of the pure substances
        from self.reference_temperature (normally 25 °C = 298.15 K) and 
        self.reference_pressure (normally 1 bar = 100000 Pa) to 
        temperature and pressure

        Parameters
        ----------
        temperature : float
            Correction temperature for the formation enthalpies. [K]
        pressure : float
            Correction pressure for the formation enthalpies. [Pa]

        """
        raise NotImplementedError()

    # ==================================================================
    # Mixtures common methods
    # ==================================================================
    def mol_fracations(self, moles: list[float]):
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
