from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt

def r_M(alpha):
  A = 2**(alpha+1)
  r = (((4*A**2 - 5*A + 1) * (alpha + 3) * zeta(alpha + 2) * zeta(alpha + 4)) / ((2*A - 1)**2 * (alpha + 2) * (zeta(alpha + 3))**2))
  return r

def dQda(alpha):
  dQda = 2*(r_M(0) - 1)/(alpha + 2)**2
  return dQda

def bisection(f, a, b, N):
  if f(a)*f(b) >= 0:
    print('Bisection method fails - 1')
    return None
  a_n = a
  b_n = b
  for n in range(1, N+1):
    print(n)
    m_n = (a_n + b_n)/2
    f_m_n = f(m_n)
    if f(a_n)*f_m_n < 0:
      a_n = a_n
      b_n = m_n
    elif f(b_n)*f_m_n < 0:
      a_n = m_n
      b_n = b_n
    elif f_m_n == 0:
      print('Found exact solution')
      return m_n
    else:
      print('Bisection method fails - 2')
      return None
  return (a_n + b_n)/2

#def zeta2(i, precision):
#  if i == 1:
#    print('Infinite')
#    exit()
#  z = 0
#  for n in range(precision):
#    z += (1/(n+1)**i)
#  return z

def get_alpha(r):
  a = -1
  b = 1
  i = 0
  j = 0
  while r < a:
    i += 1
    a = -2 + 10**(-i)
  while r > b:
    j += 1
    b = 10**(j)
  f = lambda alpha: r_M(alpha) - r
  N = 100
  alpha = bisection(f, a, b, N)
  return alpha

d = get_alpha(11)
print(d)

#alpha = np.linspace(-1.9, 50, 100)
#f = (2*r_M(0) + alpha) / (alpha + 2)

#plt.plot(alpha, f)
#plt.grid()
#plt.show()

