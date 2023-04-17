#ifndef GRID_H
#define GRID_H

struct grid
{
    int N;
    double *grid_points;
    double *dx;
    double *dx_midpoints;
    double grid_left;
    double grid_right;
};

#endif