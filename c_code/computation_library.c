#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

#include "grid.h"
#include "return_data.h"
#include "computation_data.h"
#include "physics.h"

double cal_k(double t, struct computation_data *data)
{
    return data->Lambda * exp(-t);
}

void compute_source(double t, double k, double *grid, double *u_dot, int N, struct computation_data *input_data)
{
    for (int i = 0; i < N; i++)
    {
        u_dot[i] = S(t, k, grid[i], input_data->data);
    }
}

void add_diffusion(double t, double k, double *grid, double *u, double *u_dot, int N, struct computation_data *input_data)
{
    double u_boundary;
    double *dx = input_data->computation_grid->dx;
    double *dx_midpoints = input_data->computation_grid->dx_midpoints;
    double *u_x = malloc((N + 1) * sizeof(double));
    double *Q_cal = malloc((N + 1) * sizeof(double));

    // compute u_x
    // handle left boundary
    u_boundary = left_boundary(input_data->computation_grid, u, input_data->data);
    u_x[0] = (u[0] - u_boundary) / dx[0];
    for (int i = 1; i < N; i++)
    {
        u_x[i] = (u[i] - u[i - 1]) / dx[i];
        // printf("%f\n", u_x[i]);
    }
    // handle right boundary
    u_boundary = right_boundary(input_data->computation_grid, u, input_data->data);
    u_x[N] = (u_boundary - u[N - 1]) / dx[N];

    for (int i = 0; i < N + 1; i++)
    {
        if (isnan(u_x[i]))
        {
            printf("\nnan at i=%d\n", i);
            exit(-1);
        }
        printf("\n i = %d, u_x = %f, dx = %f\n", i, u_x[i], dx[i]);
        Q_cal[i] = Q(t, k, u_x[i], input_data->data);
    }

    for (int i = 0; i < N; i++)
    {
        u_dot[i] += (Q_cal[i + 1] - Q_cal[i]) / dx_midpoints[i];
    }

    free(u_x);
    free(Q_cal);
}

int f_without_diffusion(double t, N_Vector u, N_Vector udot, void *input_data)
{
    struct computation_data *user_data = input_data;
    int N = user_data->computation_grid->N;
    double k = cal_k(t, user_data);
    double *grid = user_data->computation_grid->grid_points;
    double *udot_data = N_VGetArrayPointer(udot);

    compute_source(t, k, grid, udot_data, N, input_data);
    return 0;
}

int f_with_diffusion(double t, N_Vector u, N_Vector udot, void *input_data)
{
    struct computation_data *user_data = input_data;
    int N = user_data->computation_grid->N;
    double k = cal_k(t, user_data);
    double *grid = user_data->computation_grid->grid_points;
    double *udot_data = N_VGetArrayPointer(udot);
    double *u_data = N_VGetArrayPointer(u);

    compute_source(t, k, grid, udot_data, N, input_data);
    add_diffusion(t, k, grid, u_data, udot_data, N, input_data);
    return 0;
}

double compute_time_from_start(clock_t start_clock)
{
    return (double)(clock() - start_clock) / CLOCKS_PER_SEC;
}

void delete_previous_line()
{
    printf("\33[2K\r");
}

void log_line(double t, double k, double tir, double t_elapsed)
{
    printf("%.2f (%.4e)/%.2f - Runtime: %.1f seconds", t, k, tir, t_elapsed);
}

void compute(struct computation_data *data, struct return_data *return_struct)
{
    printf("Compuddin\n");
    sunindextype steps_to_save = return_struct->samples;
    double t_final = data->tir;
    double t_out;
    double t_now = 0.0;
    double dt = t_final / (double)(steps_to_save - 1);
    double reltol = data->tolerances;
    double abstol = data->tolerances;
    SUNContext sunctx;
    void *package_mem;
    int status;
    N_Vector u0, u_out;
    SUNMatrix jacobi_matrix;
    SUNLinearSolver lin_sol;
    SUNContext_Create(NULL, &sunctx);

    u_out = N_VNew_Serial(data->computation_grid->N, sunctx);
    u0 = N_VNew_Serial(data->computation_grid->N, sunctx);
    double *u0_c_pointer = N_VGetArrayPointer(u0);
    for (int i = 0; i < data->computation_grid->N; i++)
        u0_c_pointer[i] = initial_condition(data->computation_grid->grid_points[i], data->data);

    package_mem = CVodeCreate(CV_BDF, sunctx);

    if (activate_diffusion(data->data))
    {
        CVodeInit(package_mem, f_with_diffusion, t_now, u0);
    }
    else
    {
        CVodeInit(package_mem, f_without_diffusion, t_now, u0);
    }

    // set user data
    // ERKStepSetUserData(m_erk_mem, (void *) m_config.get())
    CVodeSetUserData(package_mem, data);

    CVodeSStolerances(package_mem, reltol, abstol);

    CVodeSetMaxNumSteps(package_mem, 1000);

    jacobi_matrix = SUNBandMatrix(data->computation_grid->N, 1, 1, sunctx);

    lin_sol = SUNLinSol_Band(u0, jacobi_matrix, sunctx);

    CVodeSetLinearSolver(package_mem, lin_sol, jacobi_matrix);

    clock_t start_clock = clock();
    log_line(t_now, cal_k(t_now, data), data->tir, compute_time_from_start(start_clock));

    // sovle shit
    for (int i = 1; i < steps_to_save; i++)
    {
        t_out = dt * i;
        do
        {
            status = CVode(package_mem, t_out, u_out, &t_now, CV_NORMAL);
            // printf("t_out = %f, t_now = %f, status = %d\n", t_out, t_now, status);
            delete_previous_line();
            log_line(t_now, cal_k(t_now, data), data->tir, compute_time_from_start(start_clock));
        } while (status == CV_TOO_MUCH_WORK);
        if (status != CV_SUCCESS)
        {
            delete_previous_line();
            printf("Error: something went wrong! CVODE error code %d\n", status);
            exit(-1);
        }
        // save after each step
        double *u_output_pointer = N_VGetArrayPointer(u_out);
        double left_point = left_boundary(data->computation_grid, u_output_pointer, data->data);
        double right_point = right_boundary(data->computation_grid, u_output_pointer, data->data);
        save_step(return_struct, i, u_output_pointer, t_now, left_point, right_point);
    }

    delete_previous_line();
    printf("Computation done in %.1f seconds", compute_time_from_start(start_clock));

    // Free stuff
    N_VDestroy(u0);
    N_VDestroy(u_out);
    SUNMatDestroy(jacobi_matrix);
    SUNLinSolFree(lin_sol);
    CVodeFree(&package_mem);
    SUNContext_Free(&sunctx);
    printf("\n");
}