import numpy as np
from Substance import Substance
from thermo.chemical import Chemical
from Mix import Mix

agua = Substance.from_thermo_database('water')
ether = Substance.from_thermo_database('ether')
benzene = Substance.from_thermo_database('benzene')
hydrogen = Substance.from_thermo_database('hydrogen')

#print(ether.name)
#print(ether.volume_gas(400,101325))
mezcla2= Mix([hydrogen, ether],'liquid')
mezcla= Mix([agua, ether, benzene], 'gas')
moles= [1000, 2000, 3000]
print (f"Las entalpias de formacion son: {mezcla.h_formations}")
print(f"entalpia agua: {benzene.h_formation}")
#print(mezcla.phase)
#print(f"concentraciones en fase {mezcla.phase}: {mezcla.concentrations(moles, 273.15, 101325)}")
cp_mezcla= mezcla.mix_heat_capacity(moles, 273.15)
print(f"La capacidad calorifica de la mezcla es: {cp_mezcla}")
print(f"Las fracciones molares son: {mezcla.mol_fracations(moles)}")
print(f"Las presiones parciales a 1 atm son: {mezcla.partial_pressures(moles, 273.15, 101325)}")

presiones_parciales = [[16887.5, 33775., 50662.5]]
print (f"Concentraciones a partir de presiones parciales: {mezcla.partial_P_2_conc(presiones_parciales,273.15)}")

# Formation enthalpies test:

print(f"Las entalpias de formacion son: {mezcla2.h_formations}")
