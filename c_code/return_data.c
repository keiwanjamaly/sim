#include "grid.h"
#include <stdlib.h>
#include <stdio.h>

struct grid *create_grid(int num_points, double *points)
{
    struct grid *new_grid = malloc(sizeof(struct grid));
    new_grid->N = num_points - 2;
    new_grid->grid_points = malloc(new_grid->N * sizeof(double));
    double *grid_midpoints = malloc((new_grid->N + 1) * sizeof(double));
    new_grid->dx = malloc((new_grid->N + 1) * sizeof(double));
    new_grid->dx_midpoints = malloc(new_grid->N * sizeof(double));

    copy_internal_grid_points(new_grid, points, num_points);
    compute_dx(new_grid);
    compute_grid_midpoints(new_grid, grid_midpoints);
    compute_dx_midpoints(new_grid, grid_midpoints);

    free(grid_midpoints);
    return new_grid;
}

void copy_internal_grid_points(struct grid *new_grid, double *points, int num_points)
{
    new_grid->grid_left = points[0];
    for (int i = 1; i < (num_points - 1); i++)
    {
        new_grid->grid_points[i - 1] = points[i];
    }
    new_grid->grid_right = points[num_points - 1];
}

void compute_dx(struct grid *new_grid)
{
    new_grid->dx[0] = new_grid->grid_points[0] - new_grid->grid_left;
    for (int i = 1; i < new_grid->N; i++)
    {
        new_grid->dx[i] = new_grid->grid_points[i] - new_grid->grid_points[i - 1];
    }
    new_grid->dx[new_grid->N] = new_grid->grid_right - new_grid->grid_points[new_grid->N - 1];
}

void compute_grid_midpoints(struct grid *new_grid, double *grid_midpoints)
{
    grid_midpoints[0] = (new_grid->grid_points[0] + new_grid->grid_left) / 2;
    for (int i = 1; i < new_grid->N; i++)
    {
        grid_midpoints[i] = (new_grid->grid_points[i] + new_grid->grid_points[i - 1]) / 2;
    }
    grid_midpoints[new_grid->N] = (new_grid->grid_right + new_grid->grid_points[new_grid->N - 1]) / 2;
}

void compute_dx_midpoints(struct grid *new_grid, double *grid_midpoints)
{
    for (int i = 0; i < new_grid->N; i++)
    {
        new_grid->dx_midpoints[i] = grid_midpoints[i + 1] - grid_midpoints[i];
    }
}

void destroy_grid(struct grid *grid_to_free)
{
    free(grid_to_free->grid_points);
    free(grid_to_free->dx);
    free(grid_to_free->dx_midpoints);
    free(grid_to_free);
}
