#include <stdio.h>

struct computation_data
{
    double Lambda;
    double kir;
    double mu;
    double T;
};

void initialize(double Lambda, double kir, double mu, double T)
{
    struct computation_data data;
    data.Lambda = Lambda;
    data.kir = kir;
    data.mu = mu;
    data.T = T;
    printf("Hello World%f, %f, %f, %f\n", data.Lambda, data.kir, data.mu, data.T);
}

void compute()
{
    printf("Compuddin\n");
}

void cleanup()
{
    printf("Cleanup\n");
}