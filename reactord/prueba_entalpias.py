from substance import Substance
from Mix import Liquid_Mix

water = Substance.from_thermo_database('naphtalene')
etanol = Substance.from_thermo_database('ethanol')
mix = Liquid_Mix([water,etanol])

print(water.normal_boiling_point)
print(water.normal_melting_point)
print(water.fusion_enthalpy(water.normal_melting_point))

print(mix.formation_enthalpies[0])
print(mix.formation_enthalpies[1])