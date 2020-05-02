// Применение метода квадратур для численного решения интегральных уравнений второго рода

import numpy as np
import math
from matplotlib import pyplot as plt
import random


def solution(x):
    #return x + np.exp(-x)                                           #1
    return 1.0                                                      #2 
    #return (25.0 + 27.0*np.cos(2*x))/(2*160*math.pi)                #3
    #return (x + ((6*x - 2 - lmbd)*lmbd)/(6 + lmbd*(lmbd - 3)))/6    #4
    #return np.cos (2*x)                                             #5
    #return 3*x*(lmbd*(2 - 3*x) + 6)/(18 + lmbd*(lmbd - 18))         #6
    #return 1 + 4*x/9                                                #7
    #return x                                                        #8
    #return 1.0                                                      #9
    #return 1.0 + 2*math.pi*lmbd*(np.cos (x)**2)/(2 - math.pi*lmbd)  #10
    #return x
    
def right_f(x):
    #return np.exp(-x)                           #1
    #return 1.0 + (np.cos(x/2) - 1)/x            #2
    #return (5.0 + 3.0*np.cos(2*x))/(16*math.pi) #3
    #return x/6                                  #4
    #return np.cos (2*x)                         #5
    #return x                                    #6
    #return 1.0 + 0.0*x                          #7
    #return 5.0*x/6.0                            #8
    #return 1.0 - x*(np.exp(x) - np.exp(-x))     #9 
    #return 1.0 + 0.0*x                          #10
    #return x + 0.5
    return 2
    
def K(x, t):
    #return x*np.exp(t)/2                                                     #1
    #return np.sin(x*t)                                                       #2
    #return -1/(4*math.pi*((np.sin(x + t)/2)**2 + 0.25*(np.cos(x + t)/2)**2)) #3
    #return lmbd*(2*x - t)                                                    #4 
    #return np.sin (x)*np.cos(t)                                              #5
    #return lmbd*x*(4*t - x)                                                  #6 
    #return x*(t**2)                                                          #7
    #return x*t/2.0                                                           #8
    #return x*x*np.exp(x*t)                                                   #9
    #return lmbd*(np.cos(x)**2) + 0.0*t                                       #10 
    return -1.0 + 0*x + 0*t

def get_approximation_values (a, b, N):

    h = (b - a)/N
    a += h/2
    
    return np.matmul (
             np.linalg.inv(np.eye (N) - np.fromfunction(lambda i, j: h*K(a + i*h, a + j*h), (N, N))),
             np.fromfunction(lambda i: right_f (a + i*h), (N, )))

def u_func (x, u, a, b):
    
    h = (b - a)/u.shape [0]
    a += h/2

    return right_f(x) + h*np.dot(u, np.fromfunction(lambda i: K(x, a + i*h), (u.shape [0], )))

def integrate (a, b, epsA, epsO, function):
    
    epsA /= b - a
    
    p = 2
    h = (b - a)*0.5
    
    result = 0.0
    while a < b:
        Ihd2 = (function (a) + function (a + h))*(h/2)
        delta = 100500.0
        while delta > epsA*h and delta > epsO*Ihd2:
            h /= 2 
            Ih = Ihd2
            Ihd21  = (function (a) + function (a + h))*(h/2)
            Ihd22 = (function (a + h) + function (a + h*2))*(h/2)
            delta = (Ih - (Ihd21 + Ihd22))/(2**p - 1)
        result += Ih
        h *= 2
        a += h
        delta = delta**p
        if delta > epsA and delta > epsO*Ihd2:
            h *= 2
        if a + h > b:
            h = b - a
    return result

def integrate_diff_K(x, y): 
    return integrate(a, b, epsA, epsO, lambda t: K(x, t) - K(y, t))

#a, b = 0.0, 1.0        #1/4/6/7/8
#a, b = 1e-6, 0.5       #2
#a, b = 1e-6, 2*math.pi #3/5
#a, b = -1.0, 1.0       #9
#a, b = 1e-6, math.pi   #10
a, b = 0.1, 1.0         #11
eps  = 1e-6
epsA = 1e-6
epsO = 1e-6
N = 20
PLOT_N = 100
lmbd = 1.0/math.pi
r1 = random.randint(1, math.trunc(10*b)) / 10


if math.fabs(integrate_diff_K(r1, b) + right_f(r1) - right_f(b)) < 2*eps:
    print('U(x) = {} - решение'.format(solution(1)))
    print("Отклонение по норме L2 от точного решения =", eps)
    print ("Отклонение по норме C[{}, {}] от точного решения =".format (math.trunc(a), b), eps)

    xList = np.linspace(a, b, PLOT_N)
    uList = np.array([solution(x) for x in xList])
    sList = np.array([solution(x) for x in xList])
    plt.plot (xList, uList, color = 'red')
    plt.plot (xList, sList, color = 'black', linestyle = 'dashed')
    plt.show()

else:    
    uNext = get_approximation_values(a, b, N)
    delta = 1.0 + eps
    while delta > eps:
    
        print ("Ошибка {0} > {1}, узлов {2}".format (delta, eps, N))
    
        N *= 2
        uPrev = uNext
        uNext = get_approximation_values(a, b, N)
        delta = np.sqrt(integrate(a, b, epsA, epsO, lambda x: (u_func(x, uNext, a, b) - u_func(x, uPrev, a, b))**2))

    print ("Функция рассчитана в {0} узлах".format (N))
    print ("Отклонение по норме в L2 от точного решения:")
    print(np.sqrt(integrate(a, b, epsA, epsO, lambda x: (u_func(x, uNext, a, b) - solution(x))**2)))
    
    xList = np.linspace(a, b, PLOT_N)
    uList = np.array([u_func (x, uNext, a, b) for x in xList])
    sList = np.array([solution (x) for x in xList])
    devList = np.abs(uList - sList)
    maxDev = devList.max ()

    print ("Отклонение по норме в C [{0}, {1}] от точного решения:".format (math.trunc(a), b))
    print (maxDev)

    plt.plot (xList, uList, color = 'green')
    plt.plot (xList, sList, color = 'black', linestyle = 'dashed')
    plt.show()
