//
// Created by Alessandro Casalino on 22/03/22.
//

#ifndef SETTINGS_H
#define SETTINGS_H

/* General parameters */
#define NEQ 2                           /* Number of equations to solve */
#define NPAR 1                          /* Number of parameters */
#define NPOINT 10000                    /* Number of points for the evolution */
#define STEPPING 0                      /* Step method - 0: linear, 1: log10 */

/* Bisection */
#define BISECTION 0                     /* Bisection target - 0: variable, 1: parameter */
#define BISECTION_TARGET 1              /* Index of the variable/parameter to use as bisection target */
#define BISECTION_MIN -1.             /* Minimum value for bisection */
#define BISECTION_MAX 5.               /* Maximum value for bisection */
#define BISECTION_PRECISION 0.01        /* Bisection precision */

/* Exit condition */
#define EXIT_CONDITION 0                /* Enable exit condition and choose the number of points to be computed after the exit condition is satisfied*/

#endif //SETTINGS_H
