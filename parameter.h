//
// Created by Alessandro Casalino on 22/03/22.
//

#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>

#include "settings.h"

// Parameters
template <typename T>
struct parameter{
    T value;
    std::string name;
};

#endif //PARAMETER_H
