#include "grid.h"
#include "unity.h"
#include "unity_fixture.h"

TEST_GROUP(GridTests);

TEST_SETUP(GridTests) {}

TEST_TEAR_DOWN(GridTests) {}

TEST(GridTests, create_and_destroy_grid) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid *test_grid = create_grid(num_points, points);

  TEST_ASSERT_NOT_NULL(test_grid);
#ifndef UNITY_EXCLUDE_DOUBLE
  TEST_ASSERT_EQUAL(num_points - 2, test_grid->N);
  TEST_ASSERT_NOT_NULL(test_grid->grid_points);
  TEST_ASSERT_NOT_NULL(test_grid->dx);
  TEST_ASSERT_NOT_NULL(test_grid->dx_midpoints);
#endif

  destroy_grid(test_grid);
}

TEST(GridTests, grid_values) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid *test_grid = create_grid(num_points, points);

  // Test internal grid points
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->grid_points[0]);
  TEST_ASSERT_EQUAL_DOUBLE(2.0, test_grid->grid_points[1]);
  TEST_ASSERT_EQUAL_DOUBLE(3.0, test_grid->grid_points[2]);

  // Test dx values
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[0]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[1]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[2]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[3]);

  // Test dx_midpoints values
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[0]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[1]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[2]);

  destroy_grid(test_grid);
}

TEST(GridTests, copy_internal_grid_points) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid test_grid;

  copy_internal_grid_points(&test_grid, num_points, points);

  TEST_ASSERT_EQUAL_DOUBLE(0.0, test_grid.grid_left);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.grid_points[0]);
  TEST_ASSERT_EQUAL_DOUBLE(2.0, test_grid.grid_points[1]);
  TEST_ASSERT_EQUAL_DOUBLE(3.0, test_grid.grid_points[2]);
  TEST_ASSERT_EQUAL_DOUBLE(4.0, test_grid.grid_right);
}

TEST(GridTests, compute_dx) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid test_grid;

  copy_internal_grid_points(&test_grid, num_points, points);
  test_grid.N = num_points - 2;
  test_grid.dx = (double *)malloc((test_grid.N + 1) * sizeof(double));

  compute_dx(&test_grid, num_points, points);

  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx[0]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx[1]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx[2]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx[3]);

  free(test_grid.dx);
}

TEST(GridTests, compute_grid_midpoints) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid test_grid;

  copy_internal_grid_points(&test_grid, num_points, points);
  test_grid.N = num_points - 2;
  double *grid_midpoints = (double *)malloc((test_grid.N + 1) * sizeof(double));

  compute_grid_midpoints(&test_grid, num_points, points, grid_midpoints);

  TEST_ASSERT_EQUAL_DOUBLE(0.5, grid_midpoints[0]);
  TEST_ASSERT_EQUAL_DOUBLE(1.5, grid_midpoints[1]);
  TEST_ASSERT_EQUAL_DOUBLE(2.5, grid_midpoints[2]);
  TEST_ASSERT_EQUAL_DOUBLE(3.5, grid_midpoints[3]);

  free(grid_midpoints);
}

TEST(GridTests, compute_dx_midpoints) {
  int num_points = 5;
  double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  struct grid test_grid;

  copy_internal_grid_points(&test_grid, num_points, points);
  test_grid.N = num_points - 2;
  double *grid_midpoints = (double *)malloc((test_grid.N + 1) * sizeof(double));
  test_grid.dx_midpoints = (double *)malloc(test_grid.N * sizeof(double));

  compute_grid_midpoints(&test_grid, num_points, points, grid_midpoints);
  compute_dx_midpoints(&test_grid, num_points, grid_midpoints);

  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx_midpoints[0]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx_midpoints[1]);
  TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid.dx_midpoints[2]);

  free(grid_midpoints);
  free(test_grid.dx_midpoints);
}

TEST_GROUP_RUNNER(GridTests) {
  RUN_TEST_CASE(GridTests, create_and_destroy_grid);
  RUN_TEST_CASE(GridTests, grid_values);
  RUN_TEST_CASE(GridTests, copy_internal_grid_points);
  RUN_TEST_CASE(GridTests, compute_dx);
  RUN_TEST_CASE(GridTests, compute_grid_midpoints);
  RUN_TEST_CASE(GridTests, compute_dx_midpoints);
}

void GridTests(void) { RUN_TEST_GROUP(GridTests); }
