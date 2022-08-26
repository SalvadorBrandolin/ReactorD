from thermo import chemical
from Substance import Substance
from Stoichiometry import Stoichiometry

Estequieometria = Stoichiometry(3,4,['hidrolisis', 'oxidacion', 'combustion'])

A=chemical.Chemical('hexane')
print(A.API)

print(A.aromatic_rings)

agua= Substance('water')

print ("\n", agua.Mw)
#print (agua.volume_s_t)

print (Estequieometria)
print (Estequieometria.coefficients)