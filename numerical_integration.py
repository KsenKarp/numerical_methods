import numpy as np
import math
import matplotlib.pyplot as plt

def f(x) :
    return (1 / (x ** 3 + x +10))

number_of_junctions = 1000

begin_x = -1
end_x = 1
h = (end_x - begin_x) / (number_of_junctions - 1)
junctions = np.linspace(begin_x, end_x, number_of_junctions)

def rect_int(f, junct, n) :
    rect_integral = 0
    for i in range (0, n - 1 ):
        rect_integral = rect_integral + ( f((junct[i + 1] + junct[i]) / 2) * (junct[i + 1] - junct[i]))
    return rect_integral
rect_integral = rect_int(f, junctions, number_of_junctions)

def trap_int(f, junct, n) :
    trap_integral = 0
    for i in range (0, n - 1 ):
        trap_integral = trap_integral + ( ((f(junct[i + 1]) + f(junct[i])) / 2) * (junct[i + 1] - junct[i]))
    return trap_integral
trap_integral = trap_int(f, junctions, number_of_junctions)

def simpson_int(f, junct, n) :
    simpson_integral = 0
    for i in range (0, n - 1 ):
        simpson_integral = simpson_integral + ((junct[i + 1] - junct[i]) / 6) * (f(junct[i]) + 4 * f((junct[i] + junct[i + 1]) / 2) + f(junct[i + 1]))
    return simpson_integral
simpson_integral = simpson_int(f, junctions, number_of_junctions)
        
print("integral calculated with rectangles method:", rect_integral)
print("integral calculated with trapeze method:", trap_integral)
print("integral calculated with Simpson's method:", simpson_integral)

def f_integral(x) :
    return np.log(abs(x+2)) / 13 - np.log(abs(x**2 - 2*x +5)) / 26 + 3 * np.arctan((x - 1) / 2) / 26

f_int_real = f_integral(1) - f_integral(-1)
print("integral real value is:", f_integral(1) - f_integral(-1))

h_log_array = []
R_log_rect = []
for j in range(number_of_junctions, number_of_junctions * 10, number_of_junctions) :
    h = (end_x - begin_x) / (j - 1)
    h_log_array.append(math.log(h))
    junctions = np.linspace(begin_x, end_x, j)
    R_log_rect.append(math.log(abs(f_int_real - rect_int(f, junctions, j))))

sp1 = plt.subplot(111)
plt.plot(h_log_array, R_log_rect, marker='.', color = 'grey')
plt.xlabel('$ln_h$', fontsize=10)
plt.grid(True)
plt.show()

R_log_trap = []
for j in range(number_of_junctions, number_of_junctions * 10, number_of_junctions) :
    #h = (end_x - begin_x) / (j - 1)
    #h_log_array.append(math.log(h))
    junctions = np.linspace(begin_x, end_x, j)
    R_log_trap.append(math.log(abs(f_int_real - trap_int(f, junctions, j))))

sp2 = plt.subplot(111)
plt.plot(h_log_array, R_log_trap, marker='.', color = 'crimson')
plt.xlabel('$ln_h$', fontsize=10)
plt.grid(True)
plt.show()
R_log_simpson = []
for j in range(number_of_junctions, number_of_junctions * 10, number_of_junctions) :
    #h = (end_x - begin_x) / (j - 1)
    #h_log_array.append(math.log(h))
    junctions = np.linspace(begin_x, end_x, j)
    R_log_simpson.append(math.log(abs(f_int_real - trap_int(f, junctions, j))))

sp3 = plt.subplot(111)
plt.plot(h_log_array, R_log_simpson, marker='.', color = 'green')
plt.xlabel('$ln_h$', fontsize=10)
plt.grid(True)
plt.show()

sp4 = plt.subplot(111)
plt.plot(h_log_array, R_log_simpson, marker='.', color = 'green')
plt.plot(h_log_array, R_log_rect, marker='.', color = 'grey')
plt.plot(h_log_array, R_log_trap, marker='.', color = 'crimson')
plt.xlabel('$ln_h$', fontsize=10)
plt.grid(True)
plt.show()

