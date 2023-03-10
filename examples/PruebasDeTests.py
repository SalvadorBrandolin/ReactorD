import numpy as np

import reactord as rd
from reactord.mix.abstract_mix import AbstractMix
from reactord.substance import Substance

r = 8.31446261815324  # m3⋅Pa/K/mol
temperature = np.array([300, 400, 500, 600])
methane = rd.Substance.from_thermo_database("methane")
oxygen = rd.Substance.from_thermo_database("oxygen")
hydrogen = rd.Substance.from_thermo_database("hydrogen")
hexane = rd.Substance.from_thermo_database("hexane")

mix1 = rd.mix.IdealGas(A=methane)
mix2 = rd.mix.IdealGas(B=oxygen)
mix3 = rd.mix.IdealGas(C=hydrogen)

# Se evalua el metodo formation_enthalpies_correction de la clase
# IdealGas (AbstractMix)


cpdt = methane.heat_capacity_gas_dt_integral(298.15, 500, 101325)
hf0 = methane.formation_enthalpy_ig
print(mix1.formation_enthalpies_correction(500, 101325))
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
mix4 = rd.mix.IdealSolution(A=hexane)
mix5 = rd.mix.IdealSolution(B=lauric_acid)
print("Ideal solutions testing:")
print(lauric_acid._sublimation_enthalpy(458.65))
print(lauric_acid._vaporization_enthalpy(460))
print("Caso con PF<298.15 -Hexano- PF= ", hexane.normal_melting_point)
print(
    "CpDt del compuesto: ",
    hexane.heat_capacity_liquid_dt_integral(298.15, 500, 101325),
)
print(
    "CpDt de la mezcla hexano: ",
    mix4.formation_enthalpies_correction(500, 101325),
    "\n\n",
)


# Caso con Punto de fusion > 298.15 K
T_Prueba = 300
PF_ac_lau = lauric_acid.normal_melting_point
Cpdt_solido = lauric_acid.heat_capacity_solid_dt_integral(
    298.15, PF_ac_lau, 101325
)
Hfusion = lauric_acid.fusion_enthalpy(PF_ac_lau)
Cpdt_liquido = lauric_acid.heat_capacity_liquid_dt_integral(
    PF_ac_lau, T_Prueba, 101325
)
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
    mix5.formation_enthalpies_correction(T_Prueba, 101325),
    "\n\n",
)

# print("H Vaporizacion: ", sucrose._vaporization_enthalpy_t(458.65))
# print(sucrose.vaporization_enthalpy(458.65))
# print(sucrose._sublimation_enthalpy_t (458.65))
# print(sucrose._vaporization_enthalpy_t (460))
# print(sucrose.fusion_enthalpy(400))


# TESTS PARA IDEAL_GAS
print("TESTS PARA IDEAL GAS\n--------------------")
co2 = rd.Substance.from_thermo_database("co2")
ethane = rd.Substance.from_thermo_database("ethane")
chlorine = rd.Substance.from_thermo_database("chlorine")

mixture = rd.mix.IdealGas(A=co2, B=ethane, C=chlorine)

print(mixture)

compositions = np.array(
    [
        [1, 1, 1],
        [10, 15, 20],
        [100, 50, 30],
        [1000, 8000, 500],
        [10000, 500, 4000],
    ]
)
temperature = np.array([300, 400, 500, 600])
pressure = np.array([101325, 150000, 200000, 300000])

# temperature = np.array([273.15])
# pressure = np.array([101325])

for t, p in zip(temperature, pressure):
    heat_cap_correction = np.array(
        [
            co2.heat_capacity_gas_dt_integral(298.15, t, p),
            ethane.heat_capacity_gas_dt_integral(298.15, t, p),
            chlorine.heat_capacity_gas_dt_integral(298.15, t, p),
        ]
    )

    raw_heat_capacities = np.array(
        [
            co2.heat_capacity_gas(t, p),
            ethane.heat_capacity_gas(t, p),
            chlorine.heat_capacity_gas(t, p),
        ]
    )

    volumes = np.array(
        [
            co2.volume_gas(t, p),
            ethane.volume_gas(t, p),
            chlorine.volume_gas(t, p),
        ]
    )
    raw_densities = p / (r * t)

    for moles in compositions:
        raw_mol_fractions = np.divide(moles, np.sum(moles))
        raw_concentrations = np.multiply(raw_mol_fractions, raw_densities)
        method_conc = mixture.concentrations(moles, t, p)
        # print ("Raw conc, Method conc: ", raw_concentrations, method_conc)

        vol_mix = mixture.volume(moles, t, p)
        vol_vikingo = sum(moles) * r * t / p
        # print("Metodo - A lo vikingo: ", vol_mix, "-", vol_vikingo)

        # Para test de mix_heat_capacity
        raw_mix_heat_capacity = np.dot(raw_heat_capacities, raw_mol_fractions)
        method_mix_heat_capacity = mixture.mix_heat_capacity(moles, t, p)
        """print(
            "raw_mix_cp: ",
            raw_mix_heat_capacity,
            "method_cp: ",
            method_mix_heat_capacity,
        )"""

        # Para test de formation_enthalpies_set
        raw_enthalpies_set = np.array(
            [
                co2.formation_enthalpy_ig,
                ethane.formation_enthalpy_ig,
                chlorine.formation_enthalpy_ig,
            ]
        )
        method_enthalpies_set = mixture._formation_enthalpies_set()
        # print(raw_enthalpies_set, method_enthalpies_set)

        # Para test de formation_enthalpies_correction
        raw_formation_enthalpies_correction = [
            co2.heat_capacity_gas_dt_integral(298.15, t, p),
            ethane.heat_capacity_gas_dt_integral(298.15, t, p),
            chlorine.heat_capacity_gas_dt_integral(298.15, t, p),
        ]

        method_formation_enthalpies_correction = (
            mixture.formation_enthalpies_correction(t, p)
        )

    """    print(
            "r_H_cor: ",
            raw_formation_enthalpies_correction,
            "\nm_H_cor: ",
            method_formation_enthalpies_correction,
            "\n",
            raw_formation_enthalpies_correction
            == method_formation_enthalpies_correction,
            "\n",
        )
        """

# ------------------------------------------------------------------------------
# PRUEBAS PARA MIXTURE_MOLECULAR_WEIGHT, MOLAR_DENSITY Y MASS_DENSITY
# ------------------------------------------------------------------------------

hexane = rd.Substance.from_thermo_database("hexane")
toluene = rd.Substance.from_thermo_database("toluene")
butanol = rd.Substance.from_thermo_database("butanol")

mixture = rd.mix.IdealSolution(A=hexane, B=toluene, C=butanol)

compositions = np.array(
    [
        [1, 5, 8],
        [10, 15, 20],
        [100, 50, 30],
        [1000, 8000, 500],
        [10000, 500, 4000],
    ]
)
temperature = np.array([300, 400, 500, 600])
pressure = 101325


print("MOLECULAR WEIGHTS")
#print(toluene.molecular_weight)
raw_molecular_weight = [
    substance.molecular_weight for substance in [hexane, toluene, butanol]
]
print(f"raw_molecular_weight: {raw_molecular_weight}")

for moles in compositions:
    raw_total_moles = sum(moles)
    raw_mol_fraction = moles / raw_total_moles
    

    # Test of mixture_molecular_weight method

    raw_mixture_molecular_weight = np.dot(
        raw_mol_fraction, raw_molecular_weight
    )

#    print(
#        f"""raw_mix_molec_weight: {np.round(raw_mixture_molecular_weight,3)},\
#    Metodo: {np.round(mixture.mixture_molecular_weight(moles),3)}, Frac molar: {raw_mol_fraction}"""
#    )

# ------------------------------------------------------------------------------
# MOLAR_DENSITY Y MASS_DENSITY
# ------------------------------------------------------------------------------
    for t in temperature:
        raw_pure_molar_volumes = np.array(
            [
                substance.volume_liquid(t, pressure)
                for substance in mixture.substances
            ]
        )

        raw_molar_density = raw_total_moles / np.dot(raw_pure_molar_volumes,moles)
        
        print(f"Raw_molar_density: {raw_molar_density}, metodo: {mixture.molar_density(moles, t, pressure)}")

        #assert np.allclose(raw_molar_density,mixture.molar_density(moles, t, pressure))

# MASS DENSITY
        raw_mass_density = raw_molar_density * raw_mixture_molecular_weight / 1000
        print(f"raw_mass_density: {raw_mass_density}, metodo: {mixture.mass_density(moles, t, pressure)}")