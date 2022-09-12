#import kinetics
import numpy as np

ri=np.array(2, 1)
ri_T=np.transpose(ri)
st=np.array([[-1, 1, 0],[-1, 2, 3]])
print(ri, ri_T)

r_c=ri_T*st

