from Substance import Substance
from thermo import Chemical
from Mix import Mix

import numpy as np

water= Substance ('water')
ether= Substance('ether')
nitrogen= Substance('nitrogen')

mezcla= Mix([water, ether, nitrogen],'liquid')
mezclaGAS= Mix([water, ether, nitrogen],'gas')

#print(f"volumen molar del N2: {nitrogen.volume_gas(273.15, 101325)*1000}")

#concentraciones= mezcla.concentrations([1000, 2000, 1500], 273.15, 101325)
#Conc_gas= mezclaGAS.concentrations([1000, 2000, 1500], 273.15, 101325)


#print(f"concentraciones en fase liquida: {concentraciones} moles/m^3")

#print(concentraciones)

#Xi= mezcla.mol_frac([1000, 2000, 1500])
#Pp= mezcla.partial_P([1000, 2000, 1500],1*101325)
#ConcFromPp= mezcla.partial_P_2_conc ([22516.66666667, 45033.33333333, 33775.], 273.15)

#print(f"Las fracciones molares son: {Xi}")
#print(f"Las presiones parciales en Pa son: {Pp}.")
#print(f"Las concentraciones a partir de las presiones parciales son: {ConcFromPp}")
#print(f"concentraciones en fase gas: {Conc_gas} moles/m^3")

mezclaX= Mix([ether, nitrogen], 'gAs')
print(mezclaX.Set_concentrations([1000, 2000, 1500]))
#print(mezclaX.concentrations)
#print(mezclaX.substances)

print (mezclaX)