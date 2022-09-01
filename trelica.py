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
def f_kb(EA,L):
    kb = np.zeros((1,4))
    for i in range(len(EA)):
        kb = np.concatenate((kb, ((EA[i]/L[i] * np.array([[1, 0, -1, 0], 
                                                            [0, 0, 0, 0], 
                                                            [-1, 0, 1, 0], 
                                                            [0, 0, 0, 0]])))))
    kb = kb[1:].T
    return kb

# Matriz de transformação das barras, T
def f_t(n, theta):
    T = np.zeros((1,4))
    for i in range(n):
        T = np.concatenate((T, 1*(np.array([[math.cos(theta[i]*math.pi/180), math.sin(theta[i]*math.pi/180), 0, 0], 
                                            [-math.sin(theta[i]*math.pi/180), math.cos(theta[i]*math.pi/180), 0, 0], 
                                            [0, 0, math.cos(theta[i]*math.pi/180), math.sin(theta[i]*math.pi/180)], 
                                            [0, 0, -math.sin(theta[i]*math.pi/180), math.cos(theta[i]*math.pi/180)]]))))
    T = T[1:].T
    return T

# Matriz de rigidez da barra no sistema global, k
def f_k(kb: np.array, T: np.array, n: int):
    k = np.zeros((1,4))

    for i in range(n):
        _t = T[:,i*4:(i+1)*4].T
        kb_ = kb[:,i*4:(i+1)*4]
        t_ = T[:,i*4:(i+1)*4]
        k = np.concatenate((k, np.matmul(_t,np.matmul(kb_,t_))))
    k = k[1:].T
    return k



#######################
# Calculo de trelicas #
#######################

def trelica(E=None, EA=None, EI=None, L=None, 
theta=None, n=None, b=None, Fe=None, P=None):
    
    # Validacao dos inputs
    args = [E,EA,EI,L,theta,n,b,Fe,P]
    x=0
    for a in args:
        if a is None:
            x += 1
    if x>2:
        print("at least 7 input arguments required")
        return

    EA = np.full((n), EA, dtype=int)

    kb = f_kb(EA,L)
    T = f_t(n, theta)
    k = f_k(kb,T,n)
    