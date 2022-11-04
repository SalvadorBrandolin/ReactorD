import numpy as np

import reactord as rd

from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance

temperature = np.array([300, 400, 500, 600])
methane = rd.Substance.from_thermo_database("methane")
oxygen = rd.Substance.from_thermo_database("oxygen")
hydrogen = rd.Substance.from_thermo_database("hydrogen")
hexane = rd.Substance.from_thermo_database("hexane")

mix1 = rd.mix.IdealGas([methane])
mix2 = rd.mix.IdealGas([oxygen])
mix3 = rd.mix.IdealGas([hydrogen])

# Se evalua el metodo formation_enthalpies_correction de la clase
# IdealGas (AbstractMix)


cpdt = methane.heat_capacity_gas_dt_integral(298.15, 500)
hf0 = methane.formation_enthalpy_ig
print(mix1.formation_enthalpies_correction(500))
print(hf0, cpdt, "hfcorregido: ", cpdt + hf0, "\n\n")

# RESULTADO: cpdt es igual a mix.1formation_enthalpies_correction(500)


class Mezcla(AbstractMix):
    def __init__(self, substance_list: list[Substance]):
        self.substances = substance_list

    def concentrations(
        self, moles: list[float], temperature: float, pressure: float
    ):
        pass

    def volume(self):
        pass

    def mix_heat_capacity(self, moles, temperature, pressure):
        pass

    def _formation_enthalpies_set(self):
        pass

    def formation_enthalpies_correction(
        self,
        temperature: float,
        pressure: float,
    ):
        pass


# -----------------------------------------------------------------------
#   Ideal solutions testing
# -----------------------------------------------------------------------

hexane = rd.Substance.from_thermo_database(
    "hexane"
)  # Sustancia con PF < 298.15
lauric_acid = rd.Substance.from_thermo_database(
    "lauric acid"
)  # Sustancia con PF > 298.15
# calcium = rd.Substance.from_thermo_database("calcium")
mix4 = rd.mix.IdealSolution([hexane])
mix5 = rd.mix.IdealSolution([lauric_acid])
print("Ideal solutions testing:")
print(lauric_acid._sublimation_enthalpy_t(458.65))
print(lauric_acid._vaporization_enthalpy_t(460))
print("Caso con PF<298.15 -Hexano- PF= ", hexane.normal_melting_point)
print(
    "CpDt del compuesto: ",
    hexane.heat_capacity_liquid_dt_integral(298.15, 500),
)
print(
    "CpDt de la mezcla hexano: ",
    mix4.formation_enthalpies_correction(500),
    "\n\n",
)


# Caso con Punto de fusion > 298.15 K
T_Prueba = 300
PF_ac_lau = lauric_acid.normal_melting_point
Cpdt_solido = lauric_acid.heat_capacity_solid_dt_integral(298.15, PF_ac_lau)
Hfusion = lauric_acid.fusion_enthalpy(PF_ac_lau)
Cpdt_liquido = lauric_acid.heat_capacity_liquid_dt_integral(PF_ac_lau, T_Prueba)
print("Caso con PF>298.15 -acido laurico- PF= ", PF_ac_lau)
print("CpDt del solido: ", Cpdt_solido)
print("entalpia de fusion: ", Hfusion)
print("CpDt del liquido: ", Cpdt_liquido, "\n")
print(
    "Correccion a traves de la sustancia: ",
    Cpdt_solido + Hfusion + Cpdt_liquido,
)

print(
    "CpDt de la mezcla acido laurico: ",
    mix5.formation_enthalpies_correction(T_Prueba),
    "\n\n",
)

# print("H Vaporizacion: ", sucrose._vaporization_enthalpy_t(458.65))
# print(sucrose.vaporization_enthalpy(458.65))

# print(sucrose._sublimation_enthalpy_t (458.65))

# print(sucrose._vaporization_enthalpy_t (460))

# print(sucrose.fusion_enthalpy(400))
