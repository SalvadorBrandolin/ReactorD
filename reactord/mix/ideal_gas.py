import numpy as np

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance


class IdealGas(AbstractMix):
    def __init__(self, substance_list: list[Substance]):
        self.substances = substance_list
        self.formation_enthalpies = [
            substance.formation_enthalpy_ig for substance in self.substances
        ]

    def concentrations(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        density = pressure / (r * temperature)
        return np.multiply(zi, density)

    def volume(self, moles, temperature, pressure):
        total_moles = np.sum(moles)

        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = total_moles * r * temperature / pressure
        return volume

    def mix_heat_capacity(self, moles, temperature, pressure):
        zi = self.mol_fracations(moles)
        pure_cp = np.array(
            [
                substance.heat_capacity_gas(temperature)
                for substance in self.substances
            ]
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
        partial_pressures = np.multiply(zi, pressure)
        return partial_pressures
    
    def formation_enthalpies(self):
        """Method that obtains/calculates the formation enthalpy of the
        mixture's substances.
        """
        raise NotImplementedError()

    def formation_enthalpies_correction(
        self,
        reference_temperature: float,
        reference_pressure: float,
        temperature: float,
        pressure: float,
    ):
        """Method that correct the formation enthalpy of the pure substances
        from self.reference_temperature (normally 25 °C = 298.15 K) and 
        self.reference_pressure (normally 1 bar = 100000 Pa) to 
        temperature and pressure.

        Parameters
        ----------
        temperature : float
            Correction temperature for the formation enthalpies. [K]
        pressure : float
            Correction pressure for the formation enthalpies. [Pa]

        """
        raise NotImplementedError()



    def partial_p_to_conc(self, partial_pressures, temperature):
        r = 8.31446261815324  # J/mol.K
        self.partial_pressures = np.array(partial_pressures)
        conc = self.partial_pressures / (r * temperature)  # mol/m^3
        return conc
