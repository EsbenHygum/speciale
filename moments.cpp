#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>

double r_M(double alpha){
  if (alpha == 1){
    return 216*log(2)*std::riemann_zeta(3)/pow(M_PI, 4);
  }
  double A = pow(2, alpha + 1);
  double r = ((4*pow(A, 2) - 5*A + 1)*(alpha + 3)*std::riemann_zeta(alpha + 2)*std::riemann_zeta(alpha + 2)) / (pow(2*A - 1, 2)*(alpha + 2)*pow(std::riemann_zeta(alpha + 3), 2));
  return r;
}

double func(double alpha, double r){
  double f = r_M(alpha) - r;
  return f;
}

double bisection(double a, double b, int N, double x_abs_tol, double epsilon, double r){
  if (func(a, r)*func(b, r) >= 0){
    std::cout << "Bisection method fails - 1" << std::endl;
    return NULL;
  }
  double a_n = a;
  double b_n = b;
  for(int x=1; x < N; x++){
    double m_n = (a_n + b_n)/2;
    double f_m_n = func(m_n, r);
    if (f_m_n == 0){
      std::cout << "Found exact solution" << std::endl;
      return m_n;
    }
    else if (func(a_n, r)*f_m_n < 0){
      b_n = m_n;
    }
    else{
      a_n = m_n;
    }
    if (abs(a_n - b_n) < x_abs_tol){
      return (a_n + b_n)/2;
    }
    if (abs(func(a_n, r) - func(b_n, r))/f_m_n < epsilon){
      return (a_n + b_n)/2;
    }
  }
  return NULL;
}

double get_alpha(double r){
  double a = -1;
  double left = a;
  double right = a;
  int i = 0;
  int j = 0;
  while (func(left, r) < 0){
    i += 1;
    right = left;
    left = -2 + pow(10, -i);
  }
  while (func(right, r) > 0){
    j += 1;
    left = right;
    right = pow(10, j);
  }
  int N = 1000;
  double tol = 1e-14;
  double epsilon = 1e-14;
  double alpha = bisection(left, right, N, tol, epsilon, r);
  return alpha;
}

int main(){
  double d = get_alpha(1.6);
  std::cout << d << std::endl;
  return 0;
}
