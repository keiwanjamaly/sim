#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "physics.h"
#include "helpers.h"

double calculate_prefactor(int dimension)
{
    if (dimension == 0)
    {
        return 0;
    }
    else
    {
        return (2 * pow(M_PI, dimension / 2.)) / (dimension * pow(2 * M_PI, dimension) * tgamma(dimension / 2.));
    }
}

double sech(double x)
{
    return 1 / cosh(x);
}

double e_b(double k, double ux)
{
    return sqrt(pow(k, 2) + ux);
}

double e_f(double k, double x)
{
    return sqrt(pow(k, 2) + pow(x, 2));
}

double n_b(double x)
{
    return 1 / expm1(x); // expm1 = exp - 1
}

double n_f(double x)
{
    return 1 / (exp(x) + 1);
}
