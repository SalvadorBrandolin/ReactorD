from Substance import Substance
from thermo.chemical import Chemical
from Mix import Liquid_Mix, IdealGas_Mix

agua = Substance.from_thermo_database('water')
ether = Substance.from_thermo_database('ether')
benzene = Substance.from_thermo_database('benzene')
hydrogen = Substance.from_thermo_database('hydrogen')


moles= [1000, 2000, 3000]

### CON CLASES ABSTRACTAS
mezcla_liq = Liquid_Mix([agua, ether, benzene])
print("MEZCLA LIQUIDA")
print(mezcla_liq)
print("Concentraciones: ", mezcla_liq.concentrations(moles, 273.15, 101325), "\n")
print("volumen: ", mezcla_liq.volume(moles, 273.15, 101325))
print("CPs: ", mezcla_liq.mix_heat_capacity(moles, 273.15, 101325))
print("Mol fractions: ", mezcla_liq.mol_fracations(moles))
print(f"Formation enthalpies: {mezcla_liq.enthalpies_formation}")

""" mezcla_gas = IdealGas_Mix([agua, ether, benzene])
print("MEZCLA GASEOSA")
print(mezcla_gas)
print("Concentraciones: ",mezcla_gas.concentrations(moles, 273.15, 101325))
print("volumen: ", mezcla_gas.volume(moles, 273.15, 101325))
print("CPs: ", mezcla_gas.mix_heat_capacity(moles, 273.15, 101325))
print("Mol fractions: ", mezcla_gas.mol_fracations(moles))
print("P parciales: ", 
      mezcla_gas.partial_pressures(moles, 273.15, 101325))
print("P parciales a concentraciones: ", 
      mezcla_gas.partial_P_2_conc ([1000,5000,6000], 298))
print(f"Formation enthalpies: {mezcla_gas.enthalpies_formation}") """


#print(f"benzene formation enthalpy: {hydrogen.h_formation}")