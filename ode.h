//
// Created by Alessandro Casalino on 21/03/22.
//

#ifndef ODE_H
#define ODE_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "settings.h"
#include "parameter.h"

template <typename T>
class ode {

public:

    // Initial condition function
    bool init_with_function{false};
    std::function<T(const T, const std::vector<T>, const std::vector<parameter<T>>)> init_var0;

    // Differential equation (in normal form)
    std::function<T(const T, const std::vector<T>, const std::vector<parameter<T>>)> diff_eq;

    // Initial condition at t_start
    double var0{0.};

    // Variable name (for log and csv)
    std::string var_name;

};


#endif //ODE_H
