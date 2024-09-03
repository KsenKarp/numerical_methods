import numpy as np
import matplotlib.pyplot as plt

def f(x) :
    return (np.arctan(np.log(x**2 + 1) + 1))**2

number_of_junctions = 100
n = number_of_junctions + 2

begin_x = -5
end_x = 5
h = (end_x - begin_x) / (number_of_junctions - 1)
begin = begin_x - h
end = end_x + h
junctions = []
for i in range (n) :
    junctions.append(begin + i*h)

#"Разлинуем" отрезок:
#junctions = np.linspace(begin, end, n) #в (-5, 5) не попадаем, нужно разлиновать именно этот отрезок и добавить две точки!

f_r_diff_array = []
for i in range (1, n-1):
    f_r_diff_array.append((f(junctions[i + 1]) - f(junctions[i])) / (junctions[i + 1] - junctions[i]))

x_right = [junctions[i] for i in range (1, n-1)]

sp1 = plt.subplot(121)
plt.plot(x_right, f_r_diff_array, marker='.', color = 'limegreen')
plt.xlabel('$x$', fontsize=10)
plt.grid(True)


f_c_diff_array = []
for i in range (1, n-1):
    f_c_diff_array.append((f(junctions[i + 1]) - f(junctions[i - 1])) / (2 * h))

x_center = [junctions[i] for i in range (1, n-1)]

plt.plot(x_center, f_c_diff_array, marker='.', color = 'crimson')
plt.xlabel('$x$', fontsize=10)
plt.legend(('$f_r(x)$', '$f_c(x), number \ of \ junctions \ is \  100$'))
plt.grid(True)

f_c_2diff_array = []
for i in range (1, n-1):
    f_c_2diff_array.append((f(junctions[i + 1]) - 2*f(junctions[i]) +f(junctions[i-1])) / (h**2))

f_c_2diff_array = []
for i in range (1, n-1):
    f_c_2diff_array.append((f(junctions[i + 1]) - 2*f(junctions[i]) +f(junctions[i-1])) / (h**2))

sp2 = plt.subplot(122)
plt.plot(x_right, f_c_2diff_array, marker='.', color = 'darkblue')
plt.xlabel('$x$', fontsize=10)
plt.grid(True)


begin1 = begin_x - 2*h
junctions_4 = []
for i in range (n+2) :
    junctions_4.append(begin1 + i*h)

f_c_2diff_array_4a = []
for i in range (2, n):
    f_c_2diff_array_4a.append((-f(junctions_4[i - 2]) + 16 * f(junctions_4[i - 1]) - 30 * f(junctions_4[i]) + 16 * f(junctions_4[i + 1]) - f(junctions_4[i + 2])) / (12 * (h**2)))

plt.plot(x_center, f_c_2diff_array_4a, marker='.', color = 'coral')
plt.xlabel('$x$', fontsize=10)
plt.legend(('$f,,_c(x), 2nd \ order\ of\ accuracy$', '$f,,_c(x), 4th \ order\ of\ accuracy$'))
plt.grid(True)
plt.show()

"""
x = sym.symbols('x')
f = (sym.atan(sym.log(x**2 + 1) + 1))**2
f_deriv = f.diff(x)
f_1_deriv = sym.lambdify(x, f_deriv, 'sympy')

R_1_r = 0
for i in range (number_of_junctions) : 
    R_1_r = R_1_r + abs((f_1_deriv(x_right[i]) - f_r_diff_array[i]))
print("R for 1st derivative calculated by right differencies: ", R_1_r)

R_1_c = 0
for i in range (number_of_junctions) : 
    R_1_c = R_1_c + abs((f_1_deriv(x_center[i]) - f_c_diff_array[i]))
print("R for 1st derivative calculated by central differencies: ", R_1_c)

f2_deriv = f.diff(x, 2)
f_2_deriv = sym.lambdify(x, f2_deriv, 'sympy')

R_2_2accuracy = 0
for i in range (number_of_junctions) : 
    R_2_2accuracy = R_2_2accuracy + abs((f_2_deriv(x_center[i]) - f_c_2diff_array[i]))
print("R for 2st derivative calculated with the 2nd order of accuracy: ", R_2_2accuracy)

R_2_4accuracy = 0
for i in range (number_of_junctions) : 
    R_2_4accuracy = R_2_4accuracy + abs((f_2_deriv(x_center[i]) - f_c_2diff_array_4a[i]))
print("R for 2st derivative calculated with the 4th order of accuracy: ", R_2_4accuracy)
"""
