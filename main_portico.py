from portico import portico
import numpy as np

##########
# Inputs #
##########

# Comprimento das barras (m)
L = np.array([3, 5, 6])

# Número de elos 
n = 4

# Módulo de elasticidade x área
EA = 6*(10**5)

EA = np.full((L.shape[0]), EA, dtype=int)# (Não alterar)

# Módulo de elasticidade x inércia
EI = 27*(10**3)

EI = np.full((L.shape[0]), EI, dtype=int)# (Não alterar)

# Ângulos das barras (graus)
theta = np.array([-90, 36.86, -90])

# Numeração dos deslocamentos nodais
b = np.array([[1,1,1,0,0,0,1,1,1,0,0,0],
              [1,1,1,1,1,1,0,0,0,0,0,0],
              [0,0,0,1,1,1,0,0,0,1,1,1]])

#Forças externas nas barras
P = np.array([[0],[40],[0]])

# Forças externas nos nós
Fe = np.array([0,0,0,0,0,-30])

##########
# Função #
##########
[U,F,f] = portico(EA,L,theta,n,b,Fe,EI,P)
print("Deslocamentos nodais global: \n",U,"\n","Forças nodais: \n",F,"\n","Forças nodais em cada barra: \n",f)