from Substance import Substance
from Mix import Mix
from Stoichiometry import Stoichiometry
from kinetics import Kinetics
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

"""Resolucion de reactor tubular con reaccion agua -> etanol"""

water = Substance.from_thermo_database('water')
ethanol = Substance.from_thermo_database('ethanol')

mix = Mix([water, ethanol], 'liquid')

estequiometria = Stoichiometry(np.array([[-1, 1]]))

#Cinética
def r1(concentrations, temperature):
    k_reaccion = 0.005 * np.exp(-15000 / (8.314 * temperature))

    return k_reaccion * concentrations[0]

cinetica = Kinetics([r1], mix, estequiometria)

#DATOS
f_in = [10, 0] #flujos de entrada
f_out = ['var', 'var'] #flujos de salida
t_in = 800 #temperatura de entrada
t_out = 'var' #temperatura de salida
dim = np.array([0, 1]) #reactor de 1 m3
a = 1 # area transversal del reactor

def q(z):           #calor como funcion de z (adiabatico)
    return 0

def mass_balance(r_i):
    return r_i

def energy_balance(r_rates, fi, temperature, q, delta_H):
    dh = np.dot(delta_H, r_rates)
    cpm = mix.mix_heat_capacity(fi, temperature)

    dT_dV = (q - dh) / (np.sum(fi) * cpm)

    return dT_dV

def border_conditions(ya, yb):
    #Condiciones de borde. En orden estas significan lo siguiente:
    #ya[0]-10 = el valor de Fagua vale 10 en la entrada (el solver va)
    # a intentar hacer que Fa-10 sea igual a cero
    bc = [ya[0]-10, ya[1], ya[2]-800]
    return bc

z = np.linspace(0, 1, 100)

def odesystem(z,var):
    fs = var[0:-1] 
    t = var[2]
    p = np.full(100, 101325)


    rs = list(map(cinetica.kinetic_eval, np.transpose(fs), t, p))
    rs = np.transpose(rs)
    r_i = rs[0]
    r_rates = rs[1]
         
    dh = list(map(cinetica.reaction_enthalpies, t, p))
        
    heat =np.transpose(list(map(q, z)))
        
    df_dz = np.transpose(list(map(mass_balance, r_i)))
    dt_dz = np.transpose(list(map(energy_balance, r_rates, np.transpose(fs), t, heat, dh)))
    
    """
    kinetic_eval = np.vectorize(cinetica._kinetic_eval_ode)
    dh = np.vectorize(cinetica.reaction_enthalpies)
    Q = np.vectorize(q)
    
    r_i, r_rates = kinetic_eval(var)
    DH = dh(t, p)
    heat = Q(z)
    
    df_dz = np.apply_along_axis(mass_balance,0, r_i)
    dt_dz = np.apply_along_axis(energy_balance, 0, r_rates, fs, t, heat, DH)
    """

    """
    kinetic_eval = np.frompyfunc(cinetica._kinetic_eval_ode, 1, 2)
    dh = np.frompyfunc(cinetica.reaction_enthalpies, 2, 1)
    Q = np.frompyfunc(q, 1, 1)

    r_i, r_rates = kinetic_eval(var)
    DH = dh(t, p)
    heat = Q(z)
    
    df_dz = np.apply_along_axis(mass_balance,0, r_i)
    dt_dz = np.apply_along_axis(energy_balance, 0, r_rates, fs, t, heat, DH)
    """

    """    
    r_i, r_rates = np.apply_along_axis(cinetica.kinetic_eval, 0, fs, t, p)
         
    dh = np.apply_along_axis(cinetica.reaction_enthalpies,0, t, p)
    heat = np.apply_along_axis(q, 0, z)
        
    df_dz = np.apply_along_axis(mass_balance,0, r_i)
    dt_dz = np.apply_along_axis(energy_balance,1, r_rates, fs, t, heat, dh)
    """
    return np.vstack((df_dz, dt_dz))

fa_guess = np.full(100, 10)
fb_guess = np.full(100, 0)
t_guess = np.full(100, 800)

y_guess = np.vstack((fa_guess, fb_guess))
y_guess = np.vstack((y_guess, t_guess))

solution = solve_bvp(odesystem, border_conditions, z, y_guess, verbose=1)

y = solution.y
z = z

plt.plot(z,y[0])
plt.show()