
import numpy as np
import matplotlib.pyplot as plt



def F_M(M, betta):
    return np.sqrt((M**2 + 2/(Gam - 1))/(2 * Gam/(Gam - 1) * M**2 * np.sin(betta)**2 - 1) +
                (M**2 * np.cos(betta)**2)/((Gam - 1)/2 * M**2 * np.sin(betta)**2 + 1))

def F_C(M, betta, c):
    return c / (((Gam - 1) * (M * np.sin(betta))**2 + 2)**(1/2) * (2 * Gam * (M * np.sin(betta))**2 - (Gam-1))**(1/2) / 
           ((Gam + 1) * M * np.sin(betta)))
        
def F_D(M, betta, d):
    return d / (((Gam + 1) * (M * np.sin(betta))**2) / ((Gam - 1) * (M * np.sin(betta))**2 + 2))

def F_P(M, betta, p):
    return p / (2 * Gam / (Gam + 1) * (M * np.sin(betta))**2 - (Gam - 1)/(Gam + 1))

def F_betta(M, betta):
    return np.arctan(2 * 1/np.tan(betta) * ((M * np.sin(betta))**2 - 1)/(M**2 * (Gam + np.cos(2 * betta)) + 2))


# Функции для определение характеристических размеров области отрыва, нормированные на площать под отрывом
def L_SR(theta_s):
    return (2 * np.sin(alpha) / (np.sin(alpha - theta_s) * np.sin(theta_s)))**(1/2)

def L_SB(theta_s, betta_s, betta_r):
    return ((np.cos(theta_s) * np.tan(theta_s + betta_s) - np.sin(theta_s)) /
            (np.cos(betta_s) * np.tan(theta_s + betta_s) - np.sin(betta_s))) * L_SR(theta_s)

def L_RB(theta_s, betta_s, betta_r):
    return ((np.cos(theta_s) * np.tan(betta_s) - np.sin(theta_s)) /
            (np.sin(theta_s + betta_r) - np.tan(betta_s) * np.cos(theta_s + betta_r))) * L_SR(theta_s)


# Функции для определения диссипации
def Eps(M_start, C_start, D_start, P_start, M_end, C_end, D_end, P_end, betta, theta):
    return 0.5 * (D_start * (M_start * C_start * np.sin(betta))**3 - D_end * (M_end * C_end * np.sin(betta-theta))**3)

def Pi(M_start, C_start, D_start, P_start, M_end, C_end, D_end, P_end, betta, theta):
    return 0.5 * (P_start + P_end) * (M_start * C_start * np.sin(betta) - M_end * C_end * np.sin(betta - theta))

def f_(M_start, C_start, D_start, P_start, M_end, C_end, D_end, P_end, betta, theta):
    return (Eps(M_start, C_start, D_start, P_start, M_end, C_end, D_end, P_end, betta, theta) - 
            Pi(M_start, C_start, D_start, P_start, M_end, C_end, D_end, P_end, betta, theta))

def F(M_0, C_0, D_0, P_0, M_1, C_1, D_1, P_1, M_2, C_2, D_2, P_2, betta_s, theta_s, betta_r, theta_r):
    return (f_(M_0, C_0, D_0, P_0, M_1, C_1, D_1, P_1, betta_s, theta_s) * L_SB(theta_s, betta_s, betta_r) +
            f_(M_1, C_1, D_1, P_1, M_2, C_2, D_2, P_2, betta_r, theta_r) * L_RB(theta_s, betta_s, betta_r))



def calc_ij(betta_s, betta_r):
    
    def catch_warning(x):
        if x < 0:
            raise Exception()
        return x
    
#     Считаем значения для области 1
    try:
        theta_s = catch_warning(F_betta(M0, betta_s))
        M1 = catch_warning(F_M(M0, betta_s))
        C1 = catch_warning(F_C(M0, betta_s, C0))
        D1 = catch_warning(F_D(M0, betta_s, D0))
        P1 = catch_warning(F_P(M0, betta_s, P0))

    #     Считаем значения для области 2
        theta_r = catch_warning(F_betta(M1, betta_r))
        if abs(theta_s + theta_r - alpha) >= eps: return np.nan, np.nan, np.nan
        if (betta_s <= theta_s or betta_r <= theta_r): return np.nan, np.nan, np.nan
        M2 = catch_warning(F_M(M1, betta_r))
        C2 = catch_warning(F_C(M1, betta_r, C1))
        D2 = catch_warning(F_D(M1, betta_r, D1))
        P2 = catch_warning(F_P(M1, betta_r, P1))

#         if (theta_r < 0 or M2 < 0 or C2 < 0 or D2 < 0 or P2 < 0):
#             return theta_s, theta_r, np.nan
    except:
        return np.nan, np.nan, np.nan
    
#     Считаем диссипацию
    F_diss = F(M0, C0, D0, P0,
               M1, C1, D1, P1,
               M2, C2, D2, P2,
               betta_s, theta_s, betta_r, theta_r)
    
    return theta_s, theta_r, F_diss

def radian_to_dergee(a):
    return a * 180 / np.pi 



# Определение конфигурации задачи
global M0, C0, D0, P0
global Gam, alpha
global eps

Gam = 1.4
alpha = 28 * np.pi / 180
eps = np.pi / (180 * 100)

M0 = 6
D0 = 1
P0 = 1
C0 = np.sqrt(Gam * P0 / D0)


# Запуск алгоритма
delta = 0.001

# задаем сетку для переменных betta
betta_s = np.arange(0.01, np.pi / 2, delta)
betta_r = np.arange(0.01, np.pi / 2, delta)
indx = range(len(betta_s))
F_diss = np.zeros([len(betta_s), len(betta_r)])
theta_s = np.zeros([len(betta_s), len(betta_r)])
theta_r = np.zeros([len(betta_s), len(betta_r)])

for i in range(len(betta_s)):
    for j in range(len(betta_r)):
        theta_s[i][j], theta_r[i][j], F_diss[i][j] = calc_ij(betta_s[i], betta_r[j])
    if i % 200 == 0:
        print("step ", i)

# Печать не nan результатов
not_nan_indxs = []

for i in range(len(betta_s)):
    for j in range(len(betta_r)):
        if np.isnan(F_diss[i][j]) == False:
            not_nan_indxs += [[i, j]]


# Печать TOP 10 не nan результатов

not_nan_F_diss = [F_diss[not_nan_indxs[n][0]][not_nan_indxs[n][1]] for n in range(len(not_nan_indxs))]
top_10_indxs = sorted(range(len(not_nan_F_diss)), key=lambda k: not_nan_F_diss[k])[:10]

for n in top_10_indxs:
    print(n, F_diss[not_nan_indxs[n][0]][not_nan_indxs[n][1]])
    print('theta_s = {0}'.format(radian_to_dergee(theta_s[not_nan_indxs[n][0]][not_nan_indxs[n][1]])))
    print('theta_r = {0}'.format(radian_to_dergee(theta_r[not_nan_indxs[n][0]][not_nan_indxs[n][1]])))
    print('betta_s = {0}'.format(radian_to_dergee(betta_s[not_nan_indxs[n][0]])))
    print('betta_r = {0}'.format(radian_to_dergee(betta_r[not_nan_indxs[n][1]])))
    print("==================================================================================")
    
n = top_10_indxs[0]
betta_best_s = betta_s[not_nan_indxs[n][0]]
betta_best_r = betta_r[not_nan_indxs[n][1]]



# Распределение F(betta_s, betta_r)

import seaborn as sns

with sns.axes_style("white"):
    f, ax = plt.subplots(figsize=(18, 10))
    ax = sns.heatmap(F_diss, 
                     vmin=25, vmax=100) 
#                      xticklabels=radian_to_dergee(betta_s), yticklabels=radian_to_dergee(betta_r))
plt.show()

