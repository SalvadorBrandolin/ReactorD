from ReactorBase import ReactorBase
import numpy as np
from Substance import Substance
from thermo.chemical import Chemical
from Mix import Mix



""" Example
A -> B reaction 1
A -> C reaction 2

[[-1,1,0],[-1,0,1]] La matriz estequiometrica
que ingreso por teclado para estas reacciones
Aunque en realidad no hace falta para probar la clase ReactorBase"""

agua = Substance.from_thermo_database('water')
ether = Substance.from_thermo_database('ether')
benzene = Substance.from_thermo_database('benzene')
mezcla = Mix([agua, ether, benzene], 'gas')

reactor_base = ReactorBase('isothermal', 'non-isobaric', mezcla)

reactor_base._r_dims_minmax = None #[[0,1],[0,1]]
reactor_base._r_volume = None #100 # List or narray containing the upper and
                              # lower bounds of the 
                              # spatial variables domain. E.g. a two 
                              # dimensions reactor: 
                              # r_dims_minmax = [[0,1],[0,1]].
# Testing the molar fluxes in/out the reactor
reactor_base._f_in = (1,'var',3)
reactor_base._f_out = ('var', 'var', 'var')


# Testing the temperature of the fluxes in/out the reactor
#reactor_base._t_in = 'var' 
#reactor_base._t_out = 'var'



# Testing the molar loads in the reactor
#reactor_base._time_minmax = None
#reactor_base._n_initial = 'var'
#reactor_base._n_final = 'var'

# Testing the Temperature in a non-isothermal reaction
#reactor_base._t_initial = None
#reactor_base._t_final = None


deg_of_freedom = reactor_base._degree_of_fredoom_check()

"""
error_count = 0
if (reactor_base._r_dims_minmax or reactor_base._r_volume) is not None:
    for i,(fin, fout) in enumerate(zip(reactor_base._f_in, reactor_base._f_out)):
        if fin == 'var' and fout == 'var': # 'var' means it's unknown
            error_count += 1
            print(error_count)
        print('Ingresó al FOR')

else:
    print('No ingresó')
"""

