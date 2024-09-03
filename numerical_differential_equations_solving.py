#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import matplotlib.pyplot as plt


# In[2]:


a = 0
b = 1
h = 0.05
h_ = 0.1
n = round((b - a) / h) + 1  # количество элементов 
x = np.linspace(a, b, n)


# In[3]:


# Исходное уравнение можно представить в виде системы уравнений
# u'(x) = y

def V(x, u, y):
    return y

def U(x, u, y):
    return (2 + (1 - 2 * x - x**2) * np.exp(x)) / (1 + x**2) + 2 * u / (1 + x**2) - 2 * x * y / (1 + x**2)

def exact(x):
    return x * np.arctan(x) - np.exp(x)


# In[4]:


# вычисляем функцию в узлах
y_ex = []
for i in x:
    y_ex.append(exact(i))


# In[5]:


# Метод Эйлера
def euler_met(x):
    u = [-1]
    y = [-1]
    for i in range(len(x) - 1):
        q0 = U(x[i], u[i], y[i])
        k0 = V(x[i], u[i], y[i])
        u.append(u[i] + h * k0)
        y.append(y[i] + h * q0)

    return u


# In[6]:


# Метод Рунге-Кутта
def runge_met(x):
    u = [-1]
    y = [-1]

    for i in range(len(x) - 1):
        q0 = U(x[i], u[i], y[i])
        k0 = V(x[i], u[i], y[i])
        q1 = U(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        k1 = V(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        q2 = U(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        k2 = V(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        q3 = U(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)
        k3 = V(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)

        u.append(u[i] + (h / 6) * (k0 + 2 * k1 + 2 * k2 + k3))
        y.append(y[i] + (h / 6) * (q0 + 2 * q1 + 2 * q2 + q3))

    return u


# In[7]:


# Метод Адамса
#логарифм ошибки от шага
def adams_met(x):
    u = [-1]
    y = [-1]

    # РК3
    for i in range(2):
        q0 = U(x[i], u[i], y[i]) * h
        k0 = V(x[i], u[i], y[i]) * h
        q1 = U(x[i] + h / 2, u[i] + k0 / 2, y[i] + q0 / 2) * h
        k1 = V(x[i] + h / 2, u[i] + k0 / 2, y[i] + q0 / 2) * h
        q2 = U(x[i] + h, u[i] - k0 + 2 * k1, y[i] - q0 + 2 * q1) * h
        k2 = V(x[i] + h, u[i] - k0 + 2 * k1, y[i] - q0 + 2 * q1) * h

        u.append(u[i] + (1 / 6) * (k0 + 4 * k1 + k2))
        y.append(y[i] + (1 / 6) * (q0 + 4 * q1 + q2))

    # метод Адамса 3-го порядка
    for i in range(2, len(x) - 1):
        q0 = U(x[i], u[i], y[i]) * h
        k0 = V(x[i], u[i], y[i]) * h
        q1 = U(x[i-1], u[i-1], y[i-1]) * h
        k1 = V(x[i-1], u[i-1], y[i-1]) * h
        q2 = U(x[i-2], u[i-2], y[i-2]) * h
        k2 = V(x[i-2], u[i-2], y[i-2]) * h

        u.append(u[i] + (1 / 12) * (23 * k0 - 16 * k1 + 5 * k2))
        y.append(y[i] + (1 / 12) * (23 * q0 - 16 * q1 + 5 * q2))

    return u


# In[8]:


# Метод Эйлера
plt.subplot(3, 1, 1)
plt.plot(x, euler_met(x), color='darkblue', label='Euler')
plt.plot(x, exact(x), color='crimson', label='exact solution') # точное решение
plt.xlabel('x', fontsize=10)
plt.ylabel('u(x)', fontsize=10)
plt.legend(loc='upper left')
plt.grid(True)

# Метод Рунге-Кутты
plt.subplot(3, 1, 2)
plt.plot(x, runge_met(x), color='darkblue', label='Runge_k')
plt.plot(x, exact(x), color='crimson', label='exact solution')
plt.xlabel('x', fontsize=10)
plt.ylabel('u(x)', fontsize=10)
plt.legend(loc='upper left')
plt.grid(True)

# Метод Адамса
plt.subplot(3, 1, 3)
plt.plot(x, adams_met(x), color='darkblue', label='Adams')
plt.plot(x, exact(x), color='crimson', label='exact solution')
plt.xlabel('x', fontsize=10)
plt.ylabel('u(x)', fontsize=10)
plt.legend(loc='upper left')
plt.grid(True)

plt.show()


# In[9]:


# погрешность метода Рунге-Кутта (4-ый порядок точности)
fault_array = []
def ln_max_err_Runge_Kutta(x, h):
    u = [-1]
    y = [-1]
    for i in range(len(x) - 1):
        q0 = U(x[i], u[i], y[i])
        k0 = V(x[i], u[i], y[i])
        q1 = U(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        k1 = V(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        q2 = U(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        k2 = V(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        q3 = U(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)
        k3 = V(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)

        u.append(u[i] + (h / 6) * (k0 + 2 * k1 + 2 * k2 + k3))
        y.append(y[i] + (h / 6) * (q0 + 2 * q1 + 2 * q2 + q3))
    err = [abs(exact(x[j]) - u[j]) for j in range(len(x))]
    fault_array.append(np.log(max(err)))


# In[10]:


ln_h = []
for i in range(3):
    n = 1 + (b - a) / h_
    x = np.linspace(a, b, int(n))

    ln_max_err_Runge_Kutta(x, h_)
    ln_h.append(np.log(h_))

    h_ /= 10

print(abs(fault_array[0] - fault_array[1]) / abs(ln_h[0] - ln_h[1]))


# In[11]:


def runge_met_(h):
    u = [0]
    y = [1]
    n = 1 + round((b - a) / h)
    x = np.linspace(a, b, n)
    

    for i in range(len(x) - 1):
        q0 = U(x[i], u[i], y[i])
        k0 = V(x[i], u[i], y[i])
        q1 = U(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        k1 = V(x[i] + h / 2, u[i] + k0 * h / 2, y[i] + q0 * h / 2)
        q2 = U(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        k2 = V(x[i] + h / 2, u[i] + k1 * h / 2, y[i] + q1 * h / 2)
        q3 = U(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)
        k3 = V(x[i] + h, u[i] + k2 * h, y[i] + q2 * h)

        u.append(u[i] + (h / 6) * (k0 + 2 * k1 + 2 * k2 + k3))
        y.append(y[i] + (h / 6) * (q0 + 2 * q1 + 2 * q2 + q3))

    return u


u = runge_met_(0.05)
u_ = runge_met_(0.1)


# In[12]:


Runge_rule = [(1 / 15) * abs(u[j] - u_[j]) for j in range(len(u_))]


# In[13]:


n = 1 + (b - a) / 0.1
x = np.linspace(a, b, int(n))
plt.plot(x, Runge_rule, color='black')
plt.xlabel('x', fontsize=10)
plt.ylabel('Runge_rule', fontsize=10)
plt.grid(True)

plt.show()


# In[14]:


fault_array = []
def ln_max_err_Adams(x, h):
    u = [-1]
    y = [-1]

    # РК3
    for i in range(2):
        q0 = U(x[i], u[i], y[i]) * h
        k0 = V(x[i], u[i], y[i]) * h
        q1 = U(x[i] + h / 2, u[i] + k0 / 2, y[i] + q0 / 2) * h
        k1 = V(x[i] + h / 2, u[i] + k0 / 2, y[i] + q0 / 2) * h
        q2 = U(x[i] + h, u[i] - k0 + 2 * k1, y[i] - q0 + 2 * q1) * h
        k2 = V(x[i] + h, u[i] - k0 + 2 * k1, y[i] - q0 + 2 * q1) * h

        u.append(u[i] + (1 / 6) * (k0 + 4 * k1 + k2))
        y.append(y[i] + (1 / 6) * (q0 + 4 * q1 + q2))

    # метод Адамса 3-го порядка
    for i in range(2, len(x) - 1):
        q0 = U(x[i], u[i], y[i]) * h
        k0 = V(x[i], u[i], y[i]) * h
        q1 = U(x[i-1], u[i-1], y[i-1]) * h
        k1 = V(x[i-1], u[i-1], y[i-1]) * h
        q2 = U(x[i-2], u[i-2], y[i-2]) * h
        k2 = V(x[i-2], u[i-2], y[i-2]) * h

        u.append(u[i] + (1 / 12) * (23 * k0 - 16 * k1 + 5 * k2))
        y.append(y[i] + (1 / 12) * (23 * q0 - 16 * q1 + 5 * q2))
    err = [abs(exact(x[j]) - u[j]) for j in range(len(x))]
    fault_array.append(np.log(max(err)))


# In[15]:


ln_h = []
for i in range(3):
    n = 1 + (b - a) / h_
    x = np.linspace(a, b, int(n))

    ln_max_err_Adams(x, h_)
    ln_h.append(np.log(h_))

    h_ /= 2

print(abs(fault_array[0] - fault_array[1]) / abs(ln_h[0] - ln_h[1]))


# In[16]:


plt.suptitle("логарифм ошибки от шага для метода Адамса", fontsize=11)

plt.plot(ln_h, fault_array, color='black')
plt.xlabel('ln(h)', fontsize=10)
plt.ylabel('ln(max_err)', fontsize=10)
plt.grid(True)

plt.show()


# In[ ]:




