# Autor: Matheus Araujo Dantas

# ============== Inputs ===============

# EA: rigidez (kN)
# EI: rigider inercial (kNm2)
# L: comprimento das barras (m)
# theta: angulação das barras (graus)
# n: número de elos
# b: numeração dos deslocamentos
# Fe: Carragamentos externos nas barras (kN ou kNm)
# P: Carregaemtnos externos nos elos (kN ou kNm)

# ======== Variáveis do programa =========

# kb: matriz de rigidez da barra no sistema local
# k: matriz de rigidez da barra no sistema global
# K: matriz de rigidez da estrutura no sistema global
# Kaa: matriz de rigidez da estrutura reduzida
# T: matriz de transformação
# Lb: matriz de incidência da barra
# u: deslocamento nodal da barra no sistema global

# =============== Outputs =============
# U: deslocamento nodal da estrutura no sistema global (mm)
# F: forças nos nós das barras (kN)
# f: forças internas nas barras no sistema local (kN)

import numpy as np
import math

###################################
# Definicao de funcoes auxiliares #
###################################

# Matriz de rigidez das barras no sistema local, kb
def f_kb(EA,EI,L):
    kb = np.zeros((6,1))
    for i in range(len(EA)):
        kb = np.hstack((kb, np.array([[EA[i]/L[i], 0, 0, -EA[i]/L[i], 0, 0], 
                                        [0, (12*EI[i])/(L[i]**3), (6*EI[i])/(L[i]**2), 0, (-12*EI[i])/(L[i]**3), (6*EI[i])/(L[i]**2)], 
                                        [0, (6*EI[i])/(L[i]**2), (4*EI[i])/(L[i]), 0, (-6*EI[i])/(L[i]**2), (2*EI[i])/(L[i])], 
                                        [-EA[i]/L[i], 0, 0, EA[i]/L[i], 0, 0],
                                        [0, (-12*EI[i])/(L[i]**3), (-6*EI[i])/(L[i]**2), 0, (12*EI[i])/(L[i]**3), (-6*EI[i])/(L[i]**2)],
                                        [0, (6*EI[i])/(L[i]**2), (2*EI[i])/(L[i]), 0, (-6*EI[i])/(L[i]**2), (4*EI[i])/(L[i])]])))
    kb = kb[:,1:]
    return kb

# Matriz de transformação das barras, T
def f_t(theta):
    T = np.zeros((6,1))
    for i in range(len(theta)):
        T = np.hstack((T, (np.array([[math.cos(theta[i]*(math.pi/180)), math.sin(theta[i]*(math.pi/180)), 0, 0, 0, 0], 
                                    [-math.sin(theta[i]*(math.pi/180)), math.cos(theta[i]*(math.pi/180)), 0, 0, 0, 0], 
                                    [0, 0, 1, 0, 0, 0],
                                    [0, 0, 0, math.cos(theta[i]*(math.pi/180)), math.sin(theta[i]*(math.pi/180))], 
                                    [0, 0, 0, -math.sin(theta[i]*(math.pi/180)), math.cos(theta[i]*(math.pi/180))],
                                    [0, 0, 0, 0, 0, 1]]))))
    T = T[:,1:]
    return T

# Matriz de rigidez da barra no sistema global, k
def f_k(kb: np.array, T: np.array, e: int):
    k = np.zeros((6,1))

    for i in range(e):
        _t = T[:,i*6:(i+1)*6].T
        kb_ = kb[:,i*4:(i+1)*6]
        t_ = T[:,i*6:(i+1)*6]
        k = np.hstack((k, np.matmul(_t,np.matmul(kb_,t_))))
    k = k[:,1:]
    return k

# Matriz de incidência, Lb
def f_Lb(e,n,b):
    Lb = np.zeros((n*3))
    for c in range(e):
        for i in  range(n*3):
            if b[c,i] == 1:
                l = np.zeros((n*3))
                l[i] = 1
                Lb = np.concatenate((Lb, l))
    Lb = Lb[n*2:]
    Lb = Lb.reshape((6*n,n*3))
    return Lb

# Matriz de rigidez da estrutura, K
def f_K(Lb: np.array, k: np.array, e: int):
    K = np.zeros((Lb.shape[1], Lb.shape[1]))
    for i in range(e):
        K = K + np.matmul(Lb[i*6:(i*6)+6,:].T, np.matmul(k[:,i*6:(i*6)+6], Lb[i*6:(i*6)+6,:]))
    return K

# Matriz de rigidez reduzida da estrutura, Kaa
def f_Kaa(K: np.array,Fe: np.array):
    l = len(Fe)
    Kaa = K[0:l,0:l]
    return Kaa

# Carregamento nodal equivalente
def f_f0(P,L,Fe):
    f0b = np.zeros(6,1)
    for i in range(len(P)):
        pass


# Deslocamentos nodais global, U [m]
def f_U(Fe,Kaa,n,f0d):
    U = np.concatenate((np.matmul(np.linalg.inv(Kaa), (Fe-f0d)),np.zeros((3*n-len(Fe))))).reshape((3*n,1))
    return U

# Forças nodais em cada barra no sistema local, f [N]
def f_f(e,kb,T,Lb,U):
    f = np.zeros((e,1))
    for i in range(e):
        f = np.concatenate((f, np.matmul(kb[:,i*6:(i*6)+6],np.matmul(T[:,i*6:(i*6)+6],np.matmul(Lb[i*6:(i*6)+6,:],U)))))
    f = f[e:]
    f = f.reshape((e,4))
    return f


#######################
# Calculo de pórticos #
#######################

def portico(EA: np.array, L: np.array, theta: np.array, 
n: int, b: np.array, Fe: np.array, EI: np.array, P: np.array):

    e = L.shape[0] # Número de barras
    kb = f_kb(EA,EI,L)
    T = f_t(theta)
    k = f_k(kb,T,e)
    Lb = f_Lb(e,n,b)
    K = f_K(Lb,k,e)
    Kaa = f_Kaa(K,Fe)
    f0d = f_f0(P,L,Fe)
    U = f_U(Fe,Kaa,n,f0d)
    F = np.matmul(K,U) # Forças nodais, F [kN]
    f = f_f(e,kb,T,Lb,U)
    return U, F, f