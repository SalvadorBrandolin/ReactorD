import numpy as np

import reactord as rd
temperature = np.array([300, 400, 500, 600])
methane = rd.Substance.from_thermo_database("methane")
oxygen = rd.Substance.from_thermo_database("oxygen")
hydrogen = rd.Substance.from_thermo_database("hydrogen")

mix1 = rd.mix.IdealGas([methane])
mix2 = rd.mix.IdealGas([oxygen])
mix3 = rd.mix.IdealGas([hydrogen])

print(mix1.formation_enthalpies_correction(500))

cpdt =  methane.heat_capacity_gas_dt_integral(298.15, 500) 
hf0 = methane.formation_enthalpy_ig
print (hf0, cpdt, "hfcorregido: ", cpdt+hf0)

# RESULTADO: cpdt es igual a mix.1formation_enthalpies_correction(500)

