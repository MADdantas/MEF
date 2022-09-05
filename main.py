import trelica as t
import numpy as np

##########
# Inputs #
##########

# Comprimento das barras (m)
L = np.array([4, 4, 3, 3, 5, 5])

# Número de elos 
n = 4

# Módulo de elasticidade x área
EA = 12*(10**4)

EA = np.full((L.shape[0]), EA, dtype=int)# (Não alterar)

# Ângulos das barras (graus)
theta = np.array([0, 0, 90, 90, 36.87, -36.87])

# Numeração dos deslocamentos nodais
b = np.array([[1,1,1,1,0,0,0,0],
              [0,0,0,0,1,1,1,1],
              [1,1,0,0,0,0,1,1],
              [0,0,1,1,1,1,0,0],
              [0,0,1,1,0,0,1,1],
              [1,1,0,0,1,1,0,0]])

# Forças externas
Fe = np.array([48, -48, 0, 0, 0])

##########
# Função #
##########
[U,F,f] = t.trelica(EA,L,theta,n,b,Fe)
print("Deslocamentos nodais global: \n",U,"\n","Forças nodais: \n",F,"\n","Forças nodais em cada barra: \n",f)
