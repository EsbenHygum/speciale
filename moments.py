def r_moments(alpha):
	precision = 1000
	A = 2**(alpha+1)
	r = ((2*A**2-5*A+1)*(alpha+3)*zeta(alpha+2, precision)*zeta(alpha+4, precision))/((2*A-1)**2*(alpha+2)*(zeta(alpha+3, precision))**2)
	return r

def zeta(i, precision):
	if i == 1:
		print('Infinite')
		exit()
	z = 0
	for n in range(precision):
		z += (1/(n+1)**i)
	return z

a = zeta(2, 10000)
print(a)

b = r_moments(1)
print(b)
