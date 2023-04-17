#include "grid.h"
#include <stdlib.h>
#include <stdio.h>

struct grid *initialize_grid(int N, double *points)
{
    struct grid *new_grid = malloc(sizeof(struct grid));
    new_grid->N = N - 2;
    new_grid->grid_points = malloc(new_grid->N * sizeof(double));
    double *grid_midpoints = malloc((new_grid->N + 1) * sizeof(double));
    new_grid->dx = malloc((new_grid->N + 1) * sizeof(double));
    new_grid->dx_midpoints = malloc(new_grid->N * sizeof(double));

    // copy grid points over
    new_grid->grid_left = points[0];
    for (int i = 1; i < (N - 1); i++)
    {
        new_grid->grid_points[i - 1] = points[i];
    }
    new_grid->grid_right = points[N - 1];

    // fill dx and grid_midpoints array
    new_grid->dx[0] = new_grid->grid_points[0] - new_grid->grid_left;
    grid_midpoints[0] = (new_grid->grid_points[0] + new_grid->grid_left) / 2;
    for (int i = 1; i < new_grid->N; i++)
    {
        new_grid->dx[i] = new_grid->grid_points[i] - new_grid->grid_points[i - 1];
        grid_midpoints[i] = (new_grid->grid_points[i] + new_grid->grid_points[i - 1]) / 2;
    }
    new_grid->dx[new_grid->N] = new_grid->grid_right - new_grid->grid_points[new_grid->N - 1];
    grid_midpoints[new_grid->N] = (new_grid->grid_right + new_grid->grid_points[new_grid->N - 1]) / 2;

    // fill dx_midpoints_array
    for (int i = 0; i < new_grid->N; i++)
    {
        new_grid->dx_midpoints[i] = grid_midpoints[i + 1] - grid_midpoints[i];
    }

    free(grid_midpoints);
    return new_grid;
}

void free_grid(struct grid *grid_to_be_freed)
{
    free(grid_to_be_freed->grid_points);
    free(grid_to_be_freed->dx);
    free(grid_to_be_freed->dx_midpoints);
    free(grid_to_be_freed);
}