import numpy as np
import matplotlib.pyplot as plt
#import sympy as sym
#from sympy.calculus.util import maximum

#интерполируемая функция
def f(x) :
    return np.sin(x/3 + np.e ** ((np.sin(x / 3))**2))

number_of_junctions = 10

begin_x = 0
end_x = 10
h = (end_x - begin_x) / (number_of_junctions - 1)
#interv = Interval(0, 10)

def Lagrange_poly(junc_array, f, x):
    s = 0
    for i in range(len(junc_array)):
        p = 1
        for j in range(len(junc_array)):
            if i != j:
                p = p * ((x - junc_array[j])/(junc_array[i]-junc_array[j]))
        s = s + p * f[i]
    return s

#"Разлинуем" отрезок:
x_equable_junctions = np.linspace(begin_x, end_x, number_of_junctions)
x_equable = [] #x* - серединки отрезков
for i in range (1, number_of_junctions):
    x_equable.append((x_equable_junctions[i-1] + x_equable_junctions[i])/2)

f_a = [f(i) for i in x_equable_junctions]
F_a = [f(i) for i in x_equable]

sp1 = plt.subplot(121)
plt.plot(x_equable, F_a, marker='.', color = 'limegreen')
plt.plot(x_equable, [Lagrange_poly(x_equable_junctions, f_a, i) for i in x_equable], linestyle='--',marker='*',  color='crimson')
plt.xlabel('$x$', fontsize=10)
plt.legend(('$f(x)$', '$L{}_n(x), number \ of \ junctions \ is \  40$'))
plt.grid(True)

x_chbsh_junctions = []
for i in range(1, number_of_junctions + 1):
    x_chbsh_junctions.append(0.5 * (begin_x + end_x) + 0.5 * (end_x - begin_x) * np.cos(
        (2 * i - 1) * np.pi / (2 * number_of_junctions)))

x_chbsh = []
for i in range(1, number_of_junctions):
    x_chbsh.append((x_chbsh_junctions[i-1] + x_chbsh_junctions[i]) / 2)

f_b = [f(i) for i in x_chbsh_junctions]
F_b = [f(i) for i in x_chbsh]

sp2 = plt.subplot(122)
plt.plot(x_chbsh, F_b, marker='.', color = 'limegreen')
plt.plot(x_chbsh, [Lagrange_poly(x_chbsh_junctions, f_b, i) for i in x_chbsh], linestyle='--',marker='*',  color='crimson')
plt.xlabel('$x$', fontsize=10)
plt.legend(('$f(x)$', '$L{}_n(x), number \ of \ junctions \ is \  40$'))
plt.grid(True)

plt.show()

for i in range(1, number_of_junctions):
    print("eq:", Lagrange_poly(x_equable_junctions, f_a, i) - f(i), end = '\n')
    print("chbsh:", Lagrange_poly(x_chbsh_junctions, f_b, i) - f(i), end = '\n')

#f_deriv = sym.diff(f, x, number_of_junctions)
#f_deriv_max = maximum(f_deriv, interv)
