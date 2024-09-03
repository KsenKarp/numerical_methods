import sympy as sym
from sympy import *
import numpy as np
import math
import matplotlib.pyplot as plt

eps_list = [10e-3, 10e-6, 10e-9]

x_begin = 0
x_end = 10

x = Symbol('x')
f = asin(tanh(x)) - x / 2 + 3

# Метод Ньютона
f_diff = diff(f)
def newton_method(f, begin, end, eps_list):
    results_new = []
    iterations_new = []
    i = 0
    for eps in eps_list:
        i = i + 1
        x_n_prev = end
        x_n = end - f.subs(x, end).n() / f_diff.subs(x, end)
        if (abs(x_n_prev - x_n) <= eps) :
            results_new.append(x_n)
            iterations_new.append(i)
        else : 
            while (abs(x_n_prev - x_n) >= eps) :
                i = i + 1
                x_n_prev = x_n
                res = (x_n - (f.subs(x, x_n) / f_diff.subs(x, x_n))).n()
                x_n = res
            results_new.append(x_n)
            iterations_new.append(i)
    return results_new, iterations_new

junctions = np.linspace(0, 10, 20)
f_b = [f.subs(x, i) for i in junctions]

x_tangent = []
x_b = 10
x_n_prev = x_b
x_n1 = (x_b - f.subs(x, x_b).n() / f_diff.subs(x, x_b)).n(9)
x_tangent.append(x_n1)

x_n_prev = x_n1
res = (x_n1 - (f.subs(x, x_n1) / f_diff.subs(x, x_n1))).n(9)
x_tangent.append(res)

k1 = ((f.subs(x, 10) - f.subs(x, x_tangent[0])).n(9)) / (10 - x_tangent[0])
b1 = f.subs(x, 10).n(9) - k1 * 10

k2 = ((f.subs(x, x_tangent[0]) - f.subs(x, x_tangent[1])).n(9)) / (x_tangent[0] - x_tangent[1])
b2 = f.subs(x, x_tangent[0]).n(9) - k1 * x_tangent[0]


f_tangent1 = [k1 * i + b1 for i in junctions]
f_tangent2 = [k2 * i + b2 for i in junctions]
sp1 = plt.subplot(111)
plt.plot(junctions, f_b, color = 'green')
plt.plot(junctions, f_tangent1, color = 'lightskyblue')
plt.plot(junctions, f_tangent2, color = 'blue')
plt.grid(True)
plt.show()
