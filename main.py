import trelica as t
import numpy as np

##########
# Inputs #
##########
EA = 12*(10**4)
L = np.array([4, 4, 3, 3, 5, 5])
n = L.shape[0]
EA = np.full((n), EA, dtype=int)
theta = np.array([0, 0, 90, 90, 36.87, -36.87])

##########
# Função #
##########
t.trelica(1,EA,3,L,theta,n,7)
