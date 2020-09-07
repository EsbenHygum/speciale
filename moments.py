from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt

def r_M(alpha):
  if (alpha == -1):
    return 216*np.log(2)*zeta(3)/np.pi**4
  A = 2**(alpha+1)
  r = (((4*A**2 - 5*A + 1) * (alpha + 3) * zeta(alpha + 2) * zeta(alpha + 4)) / ((2*A - 1)**2 * (alpha + 2) * (zeta(alpha + 3))**2))
  return r

def bisection(f, a, b, N):
  print(f(a), f(b))
  if f(a)*f(b) >= 0:
    print('Bisection method fails - 1')
    return None
  a_n = a
  b_n = b
  for n in range(1, N+1):
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

def get_alpha(r):
  a = -1
  left = a
  right = a
  i = 0
  j = 0
  f = lambda alpha: r_M(alpha) - r
  while r > f(left):
    i += 1
    right = left
    left = -2 + 10**(-i)
  while r < f(right):
    j += 1
    left = right
    right = 10**(j)
  N = 100
  alpha = bisection(f, left, right, N)
  return alpha

d = get_alpha(2)
print(d)

#r = 1.1
#alpha = np.linspace(-1.9, 50, 100)
#print(alpha[0])
#f = [0] * len(alpha)
#for ii in range(len(alpha)):
#  f[ii] = r_M(alpha[ii]) - r

#plt.plot(alpha, f)
#plt.grid()
#plt.show()

