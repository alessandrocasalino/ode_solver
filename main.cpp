#include <iostream>

#include "ode_solver.h"

int main() {

    std::cout << "ODE Solver - v.0" << std::endl;

    // Create the main ODE solver class
    ode_solver ode_solver(0.,1.);

    // Solve the differential equations
    ode_solver.rk4_solver();

    // Output the csv
    ode_solver.make_dat();

    return 0;
}
