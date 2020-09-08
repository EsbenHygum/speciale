from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt

def r_M(alpha):
  if (alpha == -1):
    return 216*np.log(2)*zeta(3)/np.pi**4
  A = 2**(alpha+1)
  r = (((4*A**2 - 5*A + 1) * (alpha + 3) * zeta(alpha + 2) * zeta(alpha + 4)) / ((2*A - 1)**2 * (alpha + 2) * (zeta(alpha + 3))**2))
  return r

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

def get_alpha(r):
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
  return alpha

d = get_alpha(1.1)
print(d)
