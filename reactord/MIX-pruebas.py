from Reactor_Mix import Mix

#Concentrations given initial moles and volume:
Mezcla= Mix(['water','methane'], phase='liquid', moles=[10,15], volume=10, temp=280, T_p=2)
print ("\n", Mezcla.concentrations, type(Mezcla.concentrations))
print ("calculo manual", 10/10, 15/10 , "\n")

##OKAY

#Concentrations given total moles, volume and molar fractions:
Mezcla2 = Mix (['chlorine', 'ether'], phase="gAs", molar_frac= [0.4, 0.6],total_moles=10, volume=15,
                T_p= 10)
print (Mezcla2.concentrations)
print ('Calculo manual:', 10*0.4/15, 10*0.6/15, "\n" )

#Partial pressures given total pressure and molar fractions:
Mezcla3 = Mix (['chlorine', 'ether'], phase="gAs",  molar_frac= [0.4, 0.6], volume=15, T_p=10)
print (Mezcla3.Par_p)
print('Calculo manual:', [0.4* 10, 0.6*10], "\n" )

#Concentrations given partial pressures:
Mezcla4= Mix (['chlorine', 'ether', 'eugenol'], phase="gAs",  Par_p= [2, 3, 4], volume=15, T_p=10, temp=350)
print (Mezcla4.concentrations)
print('Calculo manual:', [2/ (8.314*350), 3/ (8.314*350), 4/ (8.314*350)], "\n" )

#Se puede tambien inicializar el objeto directamente con las concentraciones
mezcla5= Mix(['chlorine', 'ether', 'eugenol'], phase="gAs", concentrations=[1,2,3.5])
print(mezcla5.concentrations, "Simplemente devuelve las concentraciones ingresadas")