//
// Created by Alessandro Casalino on 21/03/22.
//

#ifndef DIFF_EQS_H
#define DIFF_EQS_H

#include <vector>
#include <cmath>
#include "settings.h"
#include "ode.h"
#include "parameter.h"
#include "ode_solver.h"

/* Function definitions */
template <typename T>
extern double diff_eq_1(const T t, const std::vector<T> variables, const std::vector<parameter<T>> parameters){
    double v = variables[1];
    return v;
}

template <typename T>
extern double diff_eq_2(const T t, const std::vector<T> variables, const std::vector<parameter<T>> parameters){
    double omega = parameters[0].value;
    double x = variables[0];
    return - omega * omega * sin(x);
}

template <typename T>
extern double variable1_init(const T t, const std::vector<T> variables, const std::vector<parameter<T>> parameters){
    double x = variables[0];
    return x * x /2.;
}

#ifdef BISECTION
/* A bisection function can be specified to compute the correct initial conditions
 * The algorithm will try to find the best initial condition to satisfy the condition at t_end (or the exit condition
 * time, if specified) */
template <typename T>
extern double bisection_function(const T t, const std::vector<T> variables, const std::vector<parameter<T>> parameters){
    //double omega = parameters[0].value;
    double v = variables[1];
    return v - 1.;
}
#endif

#ifdef EXIT_CONDITION
/* An exit function can be specified to exit the computation before reaching t_end
 * Define EXIT_CONDITION > 0 to specify the number of additional points that will be
 * computed (and stored) after the exit condition is reached */
template <typename T>
extern bool exit_function(const T t, const std::vector<T> variables, const std::vector<parameter<T>> parameters){
    return variables[1] < 0.;
}
#endif

/* Communicate the definitions to the solver */
void ode_solver::communicate_to_solver(){

    /* Assign the function, the initial condition and the name of the variables */
    odes[0].diff_eq = diff_eq_1<double>;
    odes[0].var0 = M_PI_4;
    odes[0].var_name = "x";

    odes[1].diff_eq = diff_eq_2<double>;
    //odes[1].init_with_function = true;
    //odes[1].init_var0 = variable1_init<double>;
    odes[1].var_name = "v";
    odes[1].var0 = 1.;

    /* Assign the function, the initial condition and the name of the variables */
    parameters[0].value = 1.;
    parameters[0].name = "omega";

#ifdef EXIT_CONDITION
    /* Assign the exit condition function */
    exit_condition = exit_function<double>;
#endif

#ifdef BISECTION
    /* Assign the bisection target function */
    bisection_target = bisection_function<double>;
#endif

}

#endif //DIFF_EQS_H
