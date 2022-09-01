# Autor: Matheus Araujo Dantas

# ============== Inputs ===============

# E: módulo de elasticidade (kPa)
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
def f_kb(EA,L,n):
    kb = np.zeros((n,1))
    for i in range(len(EA)):
        kb = np.hstack((kb, ((EA[i]/L[i] * np.array([[1, 0, -1, 0], 
                                                            [0, 0, 0, 0], 
                                                            [-1, 0, 1, 0], 
                                                            [0, 0, 0, 0]])))))
    kb = kb[:,1:]
    return kb

# Matriz de transformação das barras, T
def f_t(theta,n,e):
    T = np.zeros((n,1))
    for i in range(len(theta)):
        T = np.hstack((T, (np.array([[math.cos(theta[i]*(math.pi/180)), math.sin(theta[i]*(math.pi/180)), 0, 0], 
                                          [-math.sin(theta[i]*(math.pi/180)), math.cos(theta[i]*(math.pi/180)), 0, 0], 
                                          [0, 0, math.cos(theta[i]*(math.pi/180)), math.sin(theta[i]*(math.pi/180))], 
                                          [0, 0, -math.sin(theta[i]*(math.pi/180)), math.cos(theta[i]*(math.pi/180))]]))))
    T = T[:,1:]
    return T

# Matriz de rigidez da barra no sistema global, k
def f_k(kb: np.array, T: np.array, e: int, n: int):
    k = np.zeros((n,1))

    for i in range(e):
        _t = T[:,i*n:(i+1)*n].T
        kb_ = kb[:,i*n:(i+1)*n]
        t_ = T[:,i*n:(i+1)*n]
        k = np.hstack((k, np.matmul(_t,np.matmul(kb_,t_))))
    k = k[:,1:]
    return k

# Matriz de incidência, Lb
def f_Lb(e,n,b):
    Lb = np.zeros((n*2))
    for c in range(e):
        for i in  range(n*2):
            if b[c,i] == 1:
                l = np.zeros((n*2))
                l[i] = 1
                Lb = np.concatenate((Lb, l))
    Lb = Lb[n*2:]
    Lb = Lb.reshape((e*n,n*2))
    return Lb

# Matriz de rigidez da estrutura, K
def f_K(Lb: np.array, k: np.array, e: int, n: int):
    K = np.zeros((Lb.shape[1], Lb.shape[1]))
    for i in range(e):
        K = K + np.matmul(Lb[i*n:(i*n)+n,:].T, np.matmul(k[:,i*n:(i*n)+n], Lb[i*n:(i*n)+n,:]))
    return K

# Matriz de rigidez reduzida da estrutura, Kaa
def f_Kaa(K: np.array,Fe: np.array):
    l = len(Fe)
    Kaa = K[0:l,0:l]
    return Kaa

# Deslocamentos nodais global, U [m]
def f_U(Fe,Kaa,n):
    U = np.concatenate((np.matmul(np.linalg.inv(Kaa), Fe),np.zeros((2*n-len(Fe))))).reshape((2*n,1))
    return U

#######################
# Calculo de trelicas #
#######################

def trelica(EA: np.array, L: np.array, 
theta: np.array, n: int, b: np.array, Fe: np.array, EI, P):
    
    # Validacao dos inputs
    args = [EA,EI,L,theta,n,b,Fe,P]
    x=0
    for a in args:
        if a is None:
            x += 1
    if x>2:
        print("at least 6 input arguments required")
        return
    e = L.shape[0] # Número de barras
    kb = f_kb(EA,L,n)
    T = f_t(theta,n,e)
    k = f_k(kb,T,e,n)
    Lb = f_Lb(e,n,b)
    K = f_K(Lb,k,e,n)
    Kaa = f_Kaa(K,Fe)
    U = f_U(Fe,Kaa,n)
    print(U)