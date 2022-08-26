import numpy as np

class Mix:
    """
    """ 

    def __init__(self, component_list, moles, phase, volume=None, Temp, Presure=None):
        self.phase= phase.lower()
        self.component_list= component_list
        self.moles= moles
        self.volume= volume
        self.concentrations= np.zeros (len(self.component_list))

        #Initial concentrations in liquid phase:
        if self.phase == 'liquid' and self.volume > 0:
            for idx, component in enumerate (self.component_list):
                self.concentrations[idx]= moles[idx] / self.volume

                #ACA SE PUEDE PONER UNA EXCEPCION SI EL VOLUMEN ES NEGATIVO O CERO (OJO CON ESTO!!)
        
    
        #Initial concentrations in gas phase:
        elif self.phase == 'gas' and self.volume > 0:
            for idx, component in enumerate (self.component_list):
                self.concentrations[idx]= Pressure[idx] / self.volume


        elif self.volume <= 0:
            print ('Volume must be a positive number') 
            


            



Mezcla= Mix(['water','methane'], [1,5], 'LiQuiD', volume=0)
print (Mezcla.concentrations, type(Mezcla.concentrations))

