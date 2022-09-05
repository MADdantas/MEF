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




#######################
# Calculo de pórticos #
#######################

def portico(EA: np.array, L: np.array, theta: np.array, 
n: int, b: np.array, Fe: np.array, EI: np.array, P: np.array):
    pass