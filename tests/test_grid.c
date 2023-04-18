#define UNITY_INCLUDE_DOUBLE
#include "unity.h"
#include "unity_fixture.h"
#include "grid.h"

// void test_create_and_destroy_grid(void)
// {
// }

TEST_GROUP(GridTests);

// extern int Counter;

TEST_SETUP(GridTests)
{
    // This is run before EACH TEST
    // Counter = 0x5a5a;
}

TEST_TEAR_DOWN(GridTests)
{
}

TEST(GridTests, create_and_destroy_grid)
{
    int num_points = 5;
    double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    struct grid *test_grid = create_grid(num_points, points);

    TEST_ASSERT_NOT_NULL(test_grid);
    TEST_ASSERT_EQUAL(num_points - 2, test_grid->N);
    TEST_ASSERT_NOT_NULL(test_grid->grid_points);
    TEST_ASSERT_NOT_NULL(test_grid->dx);
    TEST_ASSERT_NOT_NULL(test_grid->dx_midpoints);

    destroy_grid(test_grid);
}

TEST(GridTests, grid_values)
{
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

// void test_grid_values(void)
// {
//     int num_points = 5;
//     double points[] = {0.0, 1.0, 2.0, 3.0, 4.0};
//     struct grid *test_grid = create_grid(num_points, points);

//     // Test internal grid points
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->grid_points[0]);
//     TEST_ASSERT_EQUAL_DOUBLE(2.0, test_grid->grid_points[1]);
//     TEST_ASSERT_EQUAL_DOUBLE(3.0, test_grid->grid_points[2]);

//     // Test dx values
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[0]);
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[1]);
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[2]);
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx[3]);

//     // Test dx_midpoints values
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[0]);
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[1]);
//     TEST_ASSERT_EQUAL_DOUBLE(1.0, test_grid->dx_midpoints[2]);

//     destroy_grid(test_grid);
// }

// TEST_GROUP(GridTests);

// TEST_GROUP_RUNNER(GridTests)
// {
//     RUN_TEST_CASE(GridTests, test_create_and_destroy_grid);
//     RUN_TEST_CASE(GridTests, test_grid_values);
//     // RUN_TEST_CASE(GridTests, copy_internal_grid_points);
//     // RUN_TEST_CASE(GridTests, compute_dx);
//     // RUN_TEST_CASE(GridTests, compute_grid_midpoints);
//     // RUN_TEST_CASE(GridTests, compute_dx_midpoints);
// }
// Add the `void` type specifier before the TEST_GROUP macro

TEST_GROUP_RUNNER(GridTests)
{
    RUN_TEST_CASE(GridTests, create_and_destroy_grid);
    RUN_TEST_CASE(GridTests, grid_values);
    // RUN_TEST_CASE(GridTests, ComputeGridMidpoints);
    // RUN_TEST_CASE(GridTests, ComputeDxMidpoints);
}

void GridTests(void)
{
    RUN_TEST_GROUP(GridTests);
}
