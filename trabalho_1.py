import trelica as t
import numpy as np

##########
# Inputs #
##########

# Comprimento das barras (m)
L = np.array([2.5, 2.5, 2, 3.2, 2, 3.2, 2, 2.5, 2.5])

# Número de elos 
n = 6

# Módulo de elasticidade x área
EA = 237.5*(10**3)

EA = np.full((L.shape[0]), EA, dtype=int)# (Não alterar)

# Ângulos das barras (graus)
theta = np.array([0, 0, 90, 38.68, 90, 38.68, 90, 0, 0])

# Numeração dos deslocamentos nodais
b = np.array([[1,1,1,1,0,0,0,0,0,0,0,0],
              [0,0,1,1,1,1,0,0,0,0,0,0],
              [1,1,0,0,0,0,0,0,0,0,1,1],
              [0,0,1,1,0,0,0,0,0,0,1,1],
              [0,0,1,1,0,0,1,1,0,0,0,0],
              [0,0,0,0,1,1,1,1,0,0,0,0],
              [0,0,0,0,1,1,0,0,1,1,0,0],
              [0,0,0,0,0,0,1,1,0,0,1,1],
              [0,0,0,0,0,0,1,1,1,1,0,0]])

# Forças externas
Fe = np.array([0, -8.9, 0, -8.9, 0, -8.9, 0, 0, 0])

##########
# Função #
##########
[U,F,f] = t.trelica(EA,L,theta,n,b,Fe)
print("Deslocamentos nodais global: \n",U,"\n","Forças nodais: \n",F,"\n","Forças nodais em cada barra: \n",f)










