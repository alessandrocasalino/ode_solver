//
// Created by Alessandro Casalino on 21/03/22.
//

#include "ode_solver.h"
#include "definitions.h"

// Overload + operator to perform addition of vectors
template <typename T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T> & rhs){
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), [](const T & a, const T & b){ return a + b; });
    return lhs;
}

// Overload * operator to perform multiplication by a scalar
template <typename T>
std::vector<T> operator*(const T a, std::vector<T> v){
    std::transform(v.begin(), v.end(), v.begin(), [a](const T & c){ return c * a; });
    return v;
}

// Overload + to perform the addition of a vector with a scalar
// NOTE: this will add the scalar to all vector components!
template <typename T>
std::vector<T> operator+(const T a, std::vector<T> v){
    std::transform(v.begin(), v.end(), v.begin(), [a](const T & c){ return c + a; });
    return v;
}
template <typename T>
std::vector<T> operator+(std::vector<T> v, const T a){
    std::transform(v.begin(), v.end(), v.begin(), [a](const T & c){ return c + a; });
    return v;
}

// Constructor to initialize rk4 and variables
ode_solver::ode_solver(double t0, double t1) {

    // Define the starting and ending time
    t_start = t0;
    t_end = t1;

    // Fill the time vector
#if STEPPING==0
    time_linear(t_vec, t_start, t_end, NPOINT);
#else
    time_log10(t_vec, t_start, t_end, NPOINT);
#endif

    // Define the parameters and odes classes
    parameters.resize(number_parameters);
    odes.resize(number_equations);
    communicate_to_solver();

    // Initialize single-step rk4 vector
    rk4_step_vec.resize(number_equations);
    // Initialize variables vector
    variables.resize(number_equations);
    for (auto i{ variables.size() }; i-- > 0; ) {
        variables[i].resize(NPOINT);
    }

    /* Initialize variable vectors with initial conditions
     * The initialization is performed with bisection if requested */
#ifdef BISECTION
    std::cout << std::endl << "Bisection method requested for initial conditions." << std::endl << std::endl;
#if BISECTION==0
    bisection_check_variables();
#endif
    bisection(BISECTION_MIN, BISECTION_MAX, BISECTION_PRECISION);
#else
    init_conditions();
#endif

}

void ode_solver::init_conditions(const bool bisection, const double bisection_value){

    /* Firstly initialize variables with fixed value */
    for (auto i{ rk4_step_vec.size() }; i-- > 0; ) {
#ifdef BISECTION
        rk4_step_vec[i] = bisection && i==bisection_variable ? bisection_value : odes[i].var0;
#else
        rk4_step_vec[i] = odes[i].var0;
#endif
    }

    /* Then check if some variables should be initialized with functions */
    for (auto i{ rk4_step_vec.size() }; i-- > 0; ) {
        if(odes[i].init_with_function) {
            rk4_step_vec[i] = odes[i].init_var0(t_start, rk4_step_vec, parameters);
        }
    }

    /* Copy the values to the variables vector */
    for (auto i{ variables.size() }; i-- > 0; ) {
        variables[i][0] = rk4_step_vec[i];
    }

}

void ode_solver::rk4_step(double & t, const double h){

    auto size{rk4_step_vec.size()};

    // Define the vectors for the four RK4 steps
    std::vector rk0{rk4_step_vec}, rk1{rk4_step_vec}, rk2{rk4_step_vec}, rk3{rk4_step_vec};

    // Perform a RK4 step
    for (auto i{ size }; i-- > 0; ){
        rk0[i] = odes[i].diff_eq(t         , rk4_step_vec               , parameters);
    }
    for (auto i{ size }; i-- > 0; ){
        rk1[i] = odes[i].diff_eq(t + .5 * h, rk4_step_vec + .5 * h * rk0, parameters);
    }
    for (auto i{ size }; i-- > 0; ){
        rk2[i] = odes[i].diff_eq(t + .5 * h, rk4_step_vec + .5 * h * rk1, parameters);
    }
    for (auto i{ size }; i-- > 0; ){
        rk3[i] = odes[i].diff_eq(t + h     , rk4_step_vec + h * rk2     , parameters);
    }

    // Update the time and the rk4 vector
    rk4_step_vec = rk4_step_vec + 1./6. * h * (rk0 + (2.*rk1) + (2.*rk2) + rk3);

}

/* This divides an interval in a logarithmic scale of basis 10, and store results in input vector v */
void ode_solver::time_log10(std::vector<double> & v, double A, double B, int points){

    double a = log10(A);
    double b = log10(B);

    double h = (b - a) / (points-1.);

    for(int j=0; j<points; j++) {
        v.push_back( pow(10.0, a + j * h) );
    }

}

void ode_solver::time_linear(std::vector<double> & v, double A, double B, int points){

    double h = (B - A) / (points-1.);

    for(int j=0; j<points; j++) {
        v.push_back(j * h);
    }

}

double ode_solver::rk4_solver(bool bisection){

    if(!bisection){
        std::cout << "Solving the differential equation with RK4." << std::endl;
    }

    exit_index = NPOINT;

    // RK4 solver
    int j;
    for(j=0; j<exit_index-1; j++){

#ifdef EXIT_CONDITION
        // Check exit condition
        if(exit_condition(t_vec[j], rk4_step_vec, parameters) && exit_index==NPOINT) {
            exit_index = std::min(j+1+EXIT_CONDITION,NPOINT);
        }
#endif

        rk4_step(t_vec[j], t_vec[j+1]-t_vec[j]);

        for (auto i{ variables.size() }; i-- > 0; ) {
            variables[i][j+1] = rk4_step_vec[i];
        }

    }

#ifdef BISECTION
    return bisection ? bisection_target(t_vec[j-1], rk4_step_vec, parameters) : 0.;
#else
    return 0.;
#endif

}

#ifdef BISECTION
#if BISECTION==0
void ode_solver::bisection_check_variables(){
    for (auto i{ odes.size() }; i-- > 0; ) {
        if(odes[i].init_with_function && i==bisection_variable) {
            std::cout << "Warning: Can not apply bisection on variable " << odes[i].var_name << " (" << i
            << ") if its initial conditions are computed from a function. " << std::endl
            << "         Switching off initial condition computation with function." << std::endl << std::endl;
            odes[i].init_with_function = false;
        }
    }
}
#endif

bool ode_solver::bisection(double min, double max, double precision){

    double C{min};
    double rk4_min, rk4_C;

    while ( std::abs((max-min)/min) > precision ) {

#if BISECTION==1
        parameters[bisection_variable].value = min;
#else
        init_conditions(true, min);
#endif
        rk4_min = rk4_solver(true);

#if BISECTION==1
        parameters[bisection_variable].value = C;
#else
        init_conditions(true, C);
#endif
        rk4_C = rk4_solver(true);

        if(rk4_min * rk4_C >= 0) {
            min = C;
        }
        else {
            max = C;
        }

        std::cout << std::setprecision(4) << std::scientific
                  << "\tmin: " << min << " max: " << max << " C: "<< C
                  << " rk4_min: " << rk4_min << " rk4_C: " << rk4_C << std::endl;

        C = (max+min)/2.;

    }

    if(max==BISECTION_MAX || min==BISECTION_MIN){
        std::cerr << "The range used for bisection seems not to include the result." << std::endl;
        return false;
    }

    C = (max+min)/2.;

    std::cout << std::endl << "\t-> Result of bisection method is: " << C << "."
              << std::endl << std::endl;

#if BISECTION==0
    init_conditions(true, C);
#endif

    return true;
}
#endif

bool ode_solver::make_dat(){

    std::cout << "Creating the .dat file." << std::endl;

    int precision{5};
    int width{11};

    // Open file
    std::string filename{"res.dat"};
    std::ofstream file(filename, std::ios_base::out);

    if(!file.is_open()){
        std::cerr << "File cannot be created." << std::endl;
        return false;
    }

    std::string time_var_name = "time";

    // Header
    file << time_var_name << std::setw(width);
    for (int i = 0; i < number_equations; i++){
        file << odes[i].var_name << std::setw(width);
    }
    file << '\n';

    // Write data
    for(int j=0; j<exit_index-1; j++){
        file << std::setprecision(precision) << std::fixed << std::scientific << t_vec[j] << " ";
        for (int i = 0; i < number_equations; i++){
            file << std::setprecision(precision) << std::fixed << std::scientific << variables[i][j] << " ";
        }
        file << '\n';
    }

    file.close();

    return true;
}