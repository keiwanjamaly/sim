#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "physics.h"

struct physics_data
{
    double Lambda;
    double h;
    double sigma_0;
    int dimension;
    int dimension_gamma;
    double A_d;
    double N;
    double one_over_N;
    double mu;
    double beta;
};

struct physics_data *
initialize_physics_data(double Lambda, double h, double sigma_0, int dimension, int dimension_gamma, double mu, double T, double N)
{
    struct physics_data *new_physics_data = malloc(sizeof(struct physics_data));
    new_physics_data->Lambda = Lambda;
    new_physics_data->h = h;
    new_physics_data->sigma_0 = sigma_0;
    new_physics_data->dimension = dimension;
    new_physics_data->dimension_gamma = dimension_gamma;
    new_physics_data->N = N;
    new_physics_data->one_over_N = 1 / N;
    new_physics_data->mu = mu;
    new_physics_data->beta = 1 / T;
    if (dimension == 0)
    {
        new_physics_data->A_d = 0;
    }
    else
    {
        new_physics_data->A_d = (2 * pow(M_PI, dimension / 2.)) / (dimension * pow(2 * M_PI, dimension) * tgamma(dimension / 2.));
    }
    return new_physics_data;
}

void free_physics_data(struct physics_data *physics_data_to_be_freed)
{
    free(physics_data_to_be_freed);
}

bool activate_diffusion(struct physics_data *data)
{
    return !isinf(data->N);
}

double initial_condition(double x, struct physics_data *data)
{
    int d = data->dimension;
    int d_gamma = data->dimension_gamma;
    double Lambda = data->Lambda;
    double h = data->h;
    double sigma_0 = data->sigma_0;
    double tmp;
    if (d == 1)
    {
        tmp = 1 / sqrt(1.0 + pow(h * sigma_0 / Lambda, 2));
        return d_gamma * (pow(h, 2) * x / M_PI) * (atanh(tmp) - tmp);
    }
    else if (d == 2)
    {
        tmp = pow(h * sigma_0, 2) + pow(Lambda, 2);
        // TODO: check if a d_gamma needs to be inserted here
        return (2 * tmp - 2 * h * sigma_0 * sqrt(tmp)) / (2 * M_PI * sqrt(tmp)) * x * pow(h, 2);
    }
}

inline double sech(double x)
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

inline double n_b(double x)
{
    return 1 / expm1(x); // expm1 = exp - 1
}

inline double n_f(double x)
{
    return 1 / (exp(x) + 1);
}

double Q(double t, double k, double ux, struct physics_data *data)
{
    int d = data->dimension;
    int d_gamma = data->dimension_gamma;
    double A_d = data->A_d;
    double h = data->h;
    double sigma_0 = data->sigma_0;
    double e = e_b(k, ux);
    double beta = data->beta;
    double one_over_N = data->one_over_N;

    double return_value = -A_d *
                          one_over_N * pow(k, d + 2) * (1 + 2 * n_b(beta * e)) / (2 * e);
    // printf("%f\t%f\n", ux, return_value);
    if (isnan(return_value))
    {
        printf("\nnan occured ! with k=%f, ux=%f\n", k, ux);
        exit(-1);
    }
    return return_value;
}

double S(double t, double k, double x, struct physics_data *data)
{
    int d = data->dimension;
    int d_gamma = data->dimension_gamma;
    double A_d = data->A_d;
    double h = data->h;
    double sigma_0 = data->sigma_0;
    double e = e_f(k, x);
    double beta = data->beta;
    double mu = data->mu;

    double plus_mu_exponent = beta * (e + mu);
    double minus_mu_exponent = beta * (e - mu);

    double n_f_plus = n_f(plus_mu_exponent);
    double n_f_minus = n_f(minus_mu_exponent);

    double tmp = n_f_plus + n_f_minus - 1 + beta * e * (pow(sech(plus_mu_exponent * 0.5), 2) + pow(sech(minus_mu_exponent * 0.5), 2));

    double return_value = pow(h * sigma_0, 2) * pow(k, d + 2) * A_d * d_gamma * (tmp) / 2 * pow(e, 3);

    return return_value;
}

double left_boundary(struct grid *, double *u, struct physics_data *data)
{
    return 0;
}
double right_boundary(struct grid *computation_grid, double *u, struct physics_data *data)
{
    int N = computation_grid->N;
    return 2 * u[N - 1] - u[N - 2];
}