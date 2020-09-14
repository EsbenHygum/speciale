#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
using namespace std;

double r_M(double alpha){
  if (alpha == 1){
    return 216*log(2)*riemann_zeta(3)/pow(M_PI, 4);
  }
  A = pow(2, alpha + 1);
  r = ((4*pow(A, 2) - 5*A + 1)*(alpha + 3)*riemann_zeta(alpha + 2)*riemann_zeta(alpha + 2)) / (pow(2*A - 1, 2)*(alpha + 2)*pow(riemann_zeta(alpha + 3), 2));
  return r;
}

double func(double alpha, double r){
  f = r_M(alpha) - r;
  return f;
}

double bisection(double a, double b, int N, double x_abs_tol, double epsilon){
  if (func(a)*func(b) >= 0){
    cout << "Bisection method fails - 1";
    return NULL;
  }
  a_n = a;
  b_n = b;
  for(int x=1, x < N, x++){
    m_n = (a_n + b_n)/2;
    f_m_n = func(m_n);
    if (f_m_n == 0){
      cout << "Found exact solution";
      return m_n;
    }
    else if (func(a_n)*f_m_n < 0){
      b_n = m_n;
    }
    else{
      a_n = m_n;
    }
    if (abs(a_n - b_n) < x_abs_tol){
      return (a_n + b_n)/2;
    }
    if (abs(func(a_n) - func(b_n))/f_m_n < epsilon){
      return (a_n + b_n)/2;
    }
  }
  return NULL
}

double get_alpha(double r){
  double a = -1;
  double left = a;
  double right = a;
  int i = 0;
  int j = 0;
  while (func(left) < 0){
    i += 1;
    right = left;
    left = -2 + pow(10, -i);
  }
  while (func(right) > 0){
    j += 1;
    left = right;
    right = pow(10, j);
  }
  int N = 1000;
  double tol = 1e-14;
  double epsilon = 1e-14;
  double alpha = bisection(left, right, N, tol, epsilon);
  return alpha;
}

int main(){
  d = get_alpha(1.6);
  cout << d;
  return 0
}
