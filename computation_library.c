#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <time.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

// define grid functions
struct grid
{
    int N;
    double *grid_points;
};

struct grid *initialize_grid(int N, double *points)
{
    struct grid *new_grid = malloc(sizeof(struct grid));
    new_grid->N = N;
    new_grid->grid_points = malloc(N * sizeof(double));

    return new_grid;
}

void free_grid(struct grid *grid_to_be_freed)
{
    free(grid_to_be_freed->grid_points);
    free(grid_to_be_freed);
}

struct computation_data
{
    double Lambda;
    double kir;
    double tir;
    double tolerances;
    struct grid *computation_grid;
    void *data;

    // initial condition
    double (*initial_condition)(double, void *);

    // Diffusion Flux (t, k, ux)
    double (*Q)(double, double, double, void *);

    // Source (t, k, x)
    double (*S)(double, double, double, void *);

    // boundary terms
    double (*left_boundary)(struct grid *, double *, void *);
    double (*right_boundary)(struct grid *, double *, void *);
};

struct computation_data *initialize_computation_data(double Lambda, double kir, struct grid *computation_grid,
                                                     void *data, double tolerances)
{
    struct computation_data *new_computation_data = malloc(sizeof(struct computation_data));
    new_computation_data->Lambda = Lambda;
    new_computation_data->kir = kir;
    new_computation_data->tir = -log(kir / Lambda);
    new_computation_data->tolerances = tolerances;
    new_computation_data->computation_grid = computation_grid;
    new_computation_data->data = data;
    new_computation_data->Q = NULL;
    new_computation_data->initial_condition = NULL;
    new_computation_data->S = NULL;
    new_computation_data->left_boundary = NULL;
    new_computation_data->right_boundary = NULL;

    return new_computation_data;
}

void free_computation_data(struct computation_data *computation_data_to_be_freed)
{
    free(computation_data_to_be_freed);
}

struct return_data
{
    sunindextype grid_size;
    sunindextype samples;
    double *grid;
    double *solution_y;
    double *solution_time;
};

struct return_data *
initialize_return_data(sunindextype samples, struct grid *computation_grid)
{
    struct return_data *new_return_data = malloc(sizeof(struct return_data));
    new_return_data->samples = samples;
    new_return_data->grid_size = computation_grid->N;
    new_return_data->grid = (double *)malloc(computation_grid->N * sizeof(double));
    new_return_data->solution_y = (double *)malloc(computation_grid->N * samples * sizeof(double));
    new_return_data->solution_time = (double *)malloc(samples * sizeof(double));
}

void save_step(struct return_data *return_data_to_be_saved_to, sunindextype index, double *y, double time)
{
    return_data_to_be_saved_to->solution_time[index] = time;
    for (int i = 0; i < return_data_to_be_saved_to->grid_size; i++)
        return_data_to_be_saved_to->solution_y[index + i] = y[i];
}

void free_return_data(struct return_data *return_data_to_be_freed)
{
    free(return_data_to_be_freed->grid);
    free(return_data_to_be_freed->solution_y);
    free(return_data_to_be_freed->solution_time);
    free(return_data_to_be_freed);
}

double cal_k(double t, struct computation_data *data)
{
    return data->Lambda * exp(-t);
}

int f(double t, N_Vector y, N_Vector ydot, void *input_data)
{
    struct computation_data *user_data = input_data;
    double k = cal_k(t, user_data);
    double *ydata = N_VGetArrayPointer(y);
    double *fdata = N_VGetArrayPointer(ydot);
    double *grid = user_data->computation_grid->grid_points;
    sunindextype N = user_data->computation_grid->N;
    void *data = user_data->data;
    double (*S)(double, double, double, void *) = user_data->S;
    double (*Q)(double, double, double, void *) = user_data->Q;

    for (int i = 0; i < N; i++)
    {
        fdata[i] = S(t, k, grid[i], data);
    }

    // check, if a source term is defined and incloude source contribution
    if (Q)
    {
        double *u_x = malloc((N + 2) * sizeof(double));
        // handle left boundary
        double uleft = user_data->left_boundary(user_data->computation_grid, ydata, data);

        // compute u_x
        for (int i = 0; i < N + 2; i++)
        {
            u_x[i] += (ydata[i + 1] - ydata[i]); // / dx;
        }

        // for (int i = 0; i < N; i++)
        // {
        //     fdata[i] += Q(t, k, ux, data);
        // }

        free(u_x);
    }

    return 0;
}

// double compute_time_from_start(clock_t start_clock)
// {
//     return (double)(clock() - start_clock) / CLOCKS_PER_SEC;
// }

// void delete_previous_line()
// {
//     printf("\33[2K\r");
// }

// void log_line(double t, double k, double tir, double t_elapsed)
// {
//     printf("%.2f (%.4e)/%.2f - Runtime: %.1f seconds", t, k, tir, t_elapsed);
// }

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

    u0 = N_VNew_Serial(data->computation_grid->N, sunctx);
    double *u0_c_pointer = N_VGetArrayPointer(u0);
    for (int i = 0; i < data->computation_grid->N; i++)
        u0_c_pointer[i] = data->initial_condition(data->computation_grid->grid_points[i], data->data);

    package_mem = CVodeCreate(CV_BDF, sunctx);

    CVodeInit(package_mem, f, t_now, u0);

    CVodeSStolerances(package_mem, reltol, abstol);

    CVodeSetMaxNumSteps(package_mem, 1000);

    jacobi_matrix = SUNBandMatrix(data->computation_grid->N, 1, 1, sunctx);

    lin_sol = SUNLinSol_Band(u0, jacobi_matrix, sunctx);

    CVodeSetLinearSolver(package_mem, lin_sol, jacobi_matrix);

    // clock_t start_clock = clock();
    // log_line(t_now, cal_k(t_now, data), data->tir, compute_time_from_start(start_clock));

    // sovle shit
    for (int i = 1; i < steps_to_save; i++)
    {
        t_out = dt * i;
        do
        {
            status = CVode(package_mem, t_out, u_out, &t_now, CV_NORMAL);
            // delete_previous_line();
            // log_line(t_now, cal_k(t_now, data), data->tir, compute_time_from_start(start_clock));
        } while (status == CV_TOO_MUCH_WORK);
        if (status != CV_SUCCESS)
        {
            printf("Error: something went wrong! CVODE error code %d", status);
        }

        // save after each step
    }

    // Free stuff
    N_VDestroy(u0);
    N_VDestroy(u_out);
    SUNMatDestroy(jacobi_matrix);
    SUNLinSolFree(lin_sol);
    CVodeFree(package_mem);
    SUNContext_Free(&sunctx);
}

void cleanup()
{
    printf("Cleanup\n");
}