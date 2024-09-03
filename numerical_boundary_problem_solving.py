#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x):
    return x * np.arctan(x) - np.exp(x)

def exact_diff(x):
    return np.arctan(x) + x / (1 + x**2) - np.exp(x)

def p(x):
    return 2 * x / (1 + x**2)

def q(x):
    return - 2 / (1 + x**2) 

def f(x):
    return (2 + (1 - 2 * x - x**2) * np.exp(x)) / (1 + x**2)


# In[2]:


#для второго порядка логарифмическую ошибку от шага, для гамма2 больше знаков
alpha1 = 1
alpha2 = 1
beta1 = 0
beta2 = 4
gamma1 = -1
gamma2 = -9.1644

# шаг
h = 0.05

x_beginning = 0
x_ending = 1

# число узлов сетки
n = 1 + (x_ending - x_beginning) / h

# узлы
x = np.linspace(x_beginning, x_ending, int(n))


# In[3]:


def a_k(x):
    return 1 / h**2 - p(x) / (2 * h)

def b_k(x):
    return -2 / h**2 + q(x)

def c_k(x):
    return 1 / h**2 + p(x) / (2 * h)


# In[4]:


b_0 = -alpha1 / h + beta1
c_0 = alpha1 / h
f_0 = gamma1
a_n = -alpha2 / h
b_n = beta2 + alpha2 / h
f_n = gamma2

A_k_ = [-c_0 / b_0]


# In[5]:


def A_k(x):
    for i in range(1, len(x) - 1):
        A_k_.append(-c_k(x[i]) / (b_k(x[i]) + a_k(x[i]) * A_k_[i - 1]))
    A_k_.append(0)

A_k(x)

B_k_ = [f_0 / b_0]

def B_k(x):
    for i in range(1, len(x) - 1):
        B_k_.append((f(x[i]) - a_k(x[i]) * B_k_[i - 1]) / (b_k(x[i]) + a_k(x[i]) * A_k_[i - 1]))
    B_k_.append((f_n - a_n * B_k_[len(x) - 2]) / (b_n + a_n * A_k_[len(x) - 2]))

B_k(x)

u = [B_k_[len(x)-1]]


# In[6]:


def run_through_method(x):
    for i in range(len(x) - 2, -1, -1):
        u.insert(0, B_k_[i] + A_k_[i] * u[0])

run_through_method(x)


# In[7]:


plt.suptitle("Приближенное решение краевой задачи для ОДУ 1 порядок точности", fontsize=11)

plt.plot(x, u, color='darkblue', label='numerical solution')
plt.plot(x, exact_solution(x), color='crimson', label='exact solution')
plt.xlabel('x', fontsize=10)
plt.ylabel('u(x)', fontsize=10)
plt.legend(loc='upper right')
plt.grid(True)

plt.show()


# In[8]:


s = 0
for i in range (len(x)):
    s = s + abs(u[i] - exact_solution(x[i]))
    
print(s)


# In[9]:


b_0 = - 2 / h**2 + 2 * beta1 / (h * alpha1) - p(0) * beta1 / alpha1 + q(0)
c_0 = 2 / h**2
f_0 = f(0) - p(0) * gamma1 / alpha1 + 2 * gamma1 / (h * alpha1)


a_n = 2 / h**2
b_n = -2 / h**2 + q(x_ending) - 2 * beta2 / (h * alpha2) - p(x_ending) * beta2 / alpha2
f_n = f(x_ending) - p(x_ending) * gamma2 / alpha2 - 2 * gamma2 / (h * alpha2)

A_k_ = [-c_0 / b_0]


# In[10]:


def A_k(x):
    for i in range(1, len(x) - 1):
        A_k_.append(-c_k(x[i]) / (b_k(x[i]) + a_k(x[i]) * A_k_[i - 1]))
    A_k_.append(0)

A_k(x)

B_k_ = [f_0 / b_0]

def B_k(x):
    for i in range(1, len(x) - 1):
        B_k_.append((f(x[i]) - a_k(x[i]) * B_k_[i - 1]) / (b_k(x[i]) + a_k(x[i]) * A_k_[i - 1]))
    B_k_.append((f_n - a_n * B_k_[len(x) - 2]) / (b_n + a_n * A_k_[len(x) - 2]))

B_k(x)

u = [B_k_[len(x)-1]]


# In[11]:


def run_through_method(x):
    for i in range(len(x) - 2, -1, -1):
        u.insert(0, B_k_[i] + A_k_[i] * u[0])

run_through_method(x)


# In[12]:


plt.suptitle("Приближенное решение краевой задачи для ОДУ 2 порядок точности", fontsize=11)

plt.plot(x, u, color='darkblue', label='numerical solution')
plt.plot(x, exact_solution(x), color='crimson', label='exact solution')
plt.xlabel('x', fontsize=10)
plt.ylabel('u(x)', fontsize=10)
plt.legend(loc='upper right')
plt.grid(True)

plt.show()


# In[13]:


s = 0
for i in range (len(x)):
    s = s + abs(u[i] - exact_solution(x[i]))
    
print(s)


# In[14]:


gamma2 = 4 * exact_solution(1) + exact_diff(1)
print(gamma2)


# In[15]:


def a_k(x, h):
    return 1 / h**2 - p(x) / (2 * h)


def b_k(x, h):
    return -2 / h**2 + q(x)


def c_k(x, h):
    return 1 / h**2 + p(x) / (2 * h)



b_0 = -2 / h**2 + 2 * beta1 / (h * alpha1) - p(0) * beta1 / alpha1 + q(0)
c_0 = 2 / h**2
f_0 = f(0) - p(0) * gamma1 / alpha1 + 2*gamma1 / (h*alpha1)


def a_n(h):
    return -alpha2 / h


def b_n(h):
    return beta2 + alpha2 / h

def f_n(h):
    return gamma2


# In[16]:


def A_k(x, h):
    A_k_ = [-c_0 / b_0]
    for i in range(1, len(x) - 1):
        A_k_.append(-c_k(x[i], h) / (b_k(x[i], h) + a_k(x[i], h) * A_k_[i - 1]))
    A_k_.append(0)
    return A_k_


def B_k(x, h):
    A_k_ = A_k(x, h)
    B_k_ = [f_0 / b_0]
    for i in range(1, len(x) - 1):
        B_k_.append((f(x[i]) - a_k(x[i], h) * B_k_[i - 1]) / (b_k(x[i], h) + a_k(x[i], h) * A_k_[i - 1]))
    B_k_.append((f_n(h) - a_n(h) * B_k_[len(x) - 2]) / (b_n(h) + a_n(h) * A_k_[len(x) - 2]))
    return B_k_


# In[17]:


fault_array = []


def run_through_method(x, h):
    A_k_ = A_k(x, h)
    B_k_ = B_k(x, h)
    u = [B_k_[len(x) - 1]]
    for i in range(len(x) - 2, -1, -1):
        u.insert(0, B_k_[i] + A_k_[i] * u[0])

    err = [abs(exact_solution(x[j]) - u[j]) for j in range(len(x))]
    fault_array.append(np.log(max(err)))


# задание интервала
x_beginning = 0
x_ending = 1

# шаг
h = 0.05

ln_h = []

for i in range(3):
    # число узлов сетки
    
    n = 1 + (x_ending - x_beginning) / h
    # узлы
    x = np.linspace(x_beginning, x_ending, int(n))

    run_through_method(x, h)

    ln_h.append(np.log(h))

    h /= 5


# In[18]:


plt.suptitle("Оценка точности метода прогонки", fontsize=11)

plt.plot(ln_h, fault_array, color='black')
plt.xlabel('ln(h)', fontsize=10)
plt.ylabel('ln(max_err)', fontsize=10)
plt.grid(True)

plt.show()


# In[ ]:





# In[ ]:




