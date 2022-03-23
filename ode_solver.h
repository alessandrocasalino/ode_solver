//
// Created by Alessandro Casalino on 21/03/22.
//

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "settings.h"
#include "ode.h"
#include "parameter.h"

class ode_solver{

public:

    // Constructor to initialize initial conditions
    ode_solver(double t0, double t1);
    void communicate_to_solver();
    void init_conditions(bool bisection = false, double bisection_value = 0.);

    // Allocate the number of equations and parameters
    std::size_t number_equations{NEQ};
    std::size_t number_parameters{NPAR};

    // Solver parameters
    double t_start;
    double t_end;

    // Vector to store parameters properties
    std::vector<parameter<double>> parameters;
    // Vector to store ODEs properties and functions
    std::vector<ode<double>> odes;

    // RK4 solvers
    void rk4_step(double &, double);
    double rk4_solver(bool bisection = false);

    // Exit condition variables and functions
    int exit_index{NPOINT};
#ifdef EXIT_CONDITION
    std::function<bool(const double, const std::vector<double>,const std::vector<parameter<double>>)> exit_condition;
#endif

    // Bisection
#ifdef BISECTION
    void bisection_check_variables();
    bool bisection(double min, double max, double precision);
    std::size_t bisection_variable{BISECTION_TARGET};
    std::function<double(const double, const std::vector<double>,const std::vector<parameter<double>>)> bisection_target;
#endif

    // dat file maker
    bool make_dat();

private:

    // Vectors used for storing variables result
    std::vector<double> t_vec;
    std::vector<std::vector<double>> variables;

    // Vector used for one RK4 step
    std::vector<double> rk4_step_vec;

    // Time stepping creators
    static void time_log10(std::vector<double> & v, double A, double B, int points);
    static void time_linear(std::vector<double> & v, double A, double B, int points);

};

#endif //ODE_SOLVER_H
