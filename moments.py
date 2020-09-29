from scipy.special import zeta
from scipy.special import gamma
import numpy as np
import math
import matplotlib.pyplot as plt

def r_M(alpha):
  if (alpha == -1):
    return 216*np.log(2)*zeta(3)/np.pi**4
  A = 2**(alpha+1)
  r = (((4*A**2 - 5*A + 1) * (alpha + 3) * zeta(alpha + 2) * zeta(alpha + 4)) / ((2*A - 1)**2 * (alpha + 2) * (zeta(alpha + 3))**2))
  return r

def x_param(alpha, M_2, M_3):
  A = 2**(alpha+1)
  x = ((2*A - 1)*M_2*zeta(alpha + 3)*gamma(alpha + 3))/(2*(A - 1)*alpha*M_3*zeta(alpha + 2)*gamma(alpha + 2))
  return x

def q_0(alpha, M_2, M_3):
  A = 2**(alpha+1)
  q0 = ((1-A**(-1))*zeta(alpha + 2)*gamma(alpha + 2)/((alpha*x_param(alpha, M_2, M_3))**(alpha + 2)*M_2))**(1/(alpha-1))
  return q0

def bisection(f, a, b, N, x_abs_tol, epsilon):
  if f(a)*f(b) >= 0:
    print('Bisection method fails - 1')
    return None
  a_n = a
  b_n = b
  for n in range(1, N+1):
    m_n = (a_n + b_n)/2
    f_m_n = f(m_n)
    if f_m_n == 0:
      print('Found exact solution')
      return m_n
    elif f(a_n)*f_m_n < 0:
      b_n = m_n
    else:
      a_n = m_n
    if np.abs(a_n - b_n) < x_abs_tol:
      return (a_n + b_n)/2
    if np.abs((f(a_n) - f(b_n))/f_m_n) < epsilon:
      return (a_n + b_n)/2
  return None

def get_alpha(r, M_2, M_3):
  a = -1
  left = a
  right = a
  i = 0
  j = 0
  f = lambda alpha: r_M(alpha) - r
  while f(left) < 0:
    i += 1
    right = left
    left = -2 + 10**(-i)
  while f(right) > 0:
    j += 1
    left = right
    right = 10**(j)
  N = 1000
  tol = 1e-14
  epsilon = 1e-14
  alpha = bisection(f, left, right, N, tol, epsilon)
  q0 = q_0(alpha, M_2, M_3)
  x = x_param(alpha, M_2, M_3)
  return alpha, q0, x

M_2 = 90
M_3 = 90
r = 500
#alpha = []
#q0 = []
#x = []
#for ii in range(0, 100):
#  r += 0.5
#  d = get_alpha(r, M_2, M_3)
#  alpha[ii] = d[0]
#  q0[ii] = d[1]
#  x[ii] = d[2]

d = get_alpha(r, M_2, M_3)
alpha = d[0]
q0 = d[1]
x = d[2]
q = np.linspace(0.1, 10, 1000)

dist = (q/q0)**(alpha - 1)*(np.exp(alpha*x*q) + 1)**(-1)

fig, ax = plt.subplots(3)
ax[0].plot(q, q**2*dist, 'r')
ax[0].set(ylabel='$q^2f(q)$')
ax[1].plot(q, q**3*dist, 'g')
ax[1].set(ylabel='$q^3f(q)$')
ax[2].plot(q, q**4*dist, 'b')
ax[2].set(ylabel='$q^4f(q)$')
plt.xlabel('Momentum q')
ax[0].set_title('$r$ = %.2f, $M_2$ = %.2f, $M_3$ = %.2f' %(r, M_2, M_3))
ax[0].grid()
ax[1].grid()
ax[2].grid()
plt.show()

