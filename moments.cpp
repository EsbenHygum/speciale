#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <exception>

class greybody_params{
  public:
    // This function gives the parameter r_M which is the ratio M_2*M_4/(M_3^2).
    double r_M(double alpha){
      if (alpha == 1){
        return 216*log(2)*std::riemann_zeta(3)/pow(M_PI, 4);
      }
      double A = pow(2, alpha + 1);
      double r = ((4*pow(A, 2) - 5*A + 1)*(alpha + 3)*std::riemann_zeta(alpha + 2)*std::riemann_zeta(alpha + 2)) / (pow(2*A - 1, 2)*(alpha + 2)*pow(std::riemann_zeta(alpha + 3), 2));
      return r;
    }

    // This function is used in the rootfinding routine bisection.
    double func(double alpha, double r){
      double f = r_M(alpha) - r;
      return f;
    }


    // Gives the grey-body parameter x.
    double x_param(double alpha, double M_2, double M_3){
      double A = pow(2, alpha + 1);
      double x = ((2*A - 1)*M_2*std::riemann_zeta(alpha + 3)*tgamma(alpha + 3))/(2*(A - 1)*alpha*M_3*std::riemann_zeta(alpha + 2)*tgamma(alpha + 2));
      return x;
    }

    // Gives the grey-body parameter q_0.
    double q0(double alpha, double M_2, double M_3){
      double A = pow(2, alpha + 1);
      double q0 = pow(((1 - 1/A)*M_2*std::riemann_zeta(alpha + 2)*tgamma(alpha + 2))/(pow(alpha*x_param(alpha, M_2, M_3), alpha + 2)*M_2), 1/(alpha - 1));
      return q0;
    }

    // Rootfinding method using the bisection algorithm. Used to solve an equation for the grey-body parameter alpha.
    double bisection(double a, double b, int N, double x_abs_tol, double epsilon, double r){
      if (func(a, r)*func(b, r) >= 0){ // Checking whether the initial guesses is on either side of zero.
        throw std::runtime_error("Bisection method fails - initial attempt not bounded correctly.");
      }
      double a_n = a;
      double b_n = b;
      for(int x=1; x < N; x++){
        double m_n = (a_n + b_n)/2; // Midway point.
        double f_m_n = func(m_n, r);
        if (f_m_n == 0){
          std::cout << "Found exact solution" << std::endl;
          return m_n;
        }
        else if (func(a_n, r)*f_m_n < 0){ 
          b_n = m_n;
        }
        else{
          a_n = m_n; // Either moving the midway point to a or b, depending on which still encapsulates zero.
        }
        if (abs(a_n - b_n) < x_abs_tol){
          return (a_n + b_n)/2;
        }
        if (abs(func(a_n, r) - func(b_n, r))/f_m_n < epsilon){
          return (a_n + b_n)/2; // Once it's close enough to zero, by some metric, return that value.
        }
      }
      throw std::runtime_error("Bisection method fails - unable to converge to a solution.");
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
      } // Searching for the left boundary, dragging the right boundary along. Guesses closer and closer to -2.
      while (func(right, r) > 0){
        j += 1;
        left = right;
        right = pow(10, j);
      } // Right boundary.
      int N = 1000; // Max iterations of the bisection routine. Should be high enough to allow for the bisection method to converge to a solution.
      double tol = 1e-14;
      double epsilon = 1e-14; // tol and epsilon is how close the bisection method comes to zero before it accepts it as a solution.
      double alpha = bisection(left, right, N, tol, epsilon, r);
    return alpha;
    }

    greybody_params(double r, double M_2, double M_3){
      double alpha = get_alpha(r);
      double x = x_param(alpha, M_2, M_3);
      double q_0 = q0(alpha, M_2, M_3);
    }
};

int main(){
  greybody_params test(1.6, 2, 2);
  std::cout << test.get_alpha(1.6) << std::endl;
  return 0;
}

