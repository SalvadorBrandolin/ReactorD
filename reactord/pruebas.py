from thermo import chemical
from Substance import Substance

A=chemical.Chemical('hexane')
print(A.API)

print(A.aromatic_rings)

agua= Substance('water')

print ("\n", agua.mw)
print (agua.volume_s_t)