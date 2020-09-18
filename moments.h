#ifndef GREYBODYPARAMS_H
#define GREYBODYPARAMS_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <exception>


class greybodyparams{
  public:
    greybodyparams(double r, double M_2, double M_3);
    const double q0_;
    const double x_;
    const double alpha_;

  private:
    // This function gives the parameter r_M which is the ratio M_2*M_4/(M_3^2).
    double r_M(double alpha);

    // This function is used in the rootfinding routine bisection.
    double func(double alpha, double r);

    // Gives the grey-body parameter x.
    double x_param(double alpha, double M_2, double M_3);

    // Gives the grey-body parameter q_0.
    double q0(double alpha, double M_2, double M_3);

    // Rootfinding method using the bisection algorithm. Used to solve an equation for the grey-body parameter alpha.
    double bisection(double a, double b, int N, double x_abs_tol, double epsilon, double r);

    double get_alpha(double r);
};

#endif
