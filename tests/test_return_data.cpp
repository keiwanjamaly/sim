#include "unity.h"
#include "unity_fixture.h"
#include "return_data.h"
#include "grid.h"

TEST_GROUP(ReturnDataTests);

TEST_SETUP(ReturnDataTests)
{
}

TEST_TEAR_DOWN(ReturnDataTests)
{
}

TEST(ReturnDataTests, create_and_destroy_grid)
{
    int samples = 3;
    int N = 4;
    double grid_points[] = {0.25, 0.5, 0.75, 1.0};
    Grid *computation_grid = create_grid(N, grid_points);

    ReturnData *result = create_return_data(samples, computation_grid);

    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(samples, result->samples);
    TEST_ASSERT_EQUAL_INT(N, result->grid_size);

    TEST_ASSERT_NOT_NULL(result->grid);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(grid_points, result->grid, N);

    TEST_ASSERT_NOT_NULL(result->solution_y);
    for (int i = 0; i < samples; i++)
    {
        TEST_ASSERT_NOT_NULL(result->solution_y[i]);
    }

    TEST_ASSERT_NOT_NULL(result->solution_time);

    // Don't forget to clean up
    destroy_return_data(result);
    destroy_grid(computation_grid);
}

TEST(ReturnDataTests, save_step_basic)
{
    int samples = 2;
    int N = 4;
    double points[] = {0.25, 0.5, 0.75, 1.0};

    Grid *computation_grid = create_grid(N, points);

    ReturnData *result = create_return_data(samples, computation_grid);

    int index = 0;
    double y[] = {1.0, 2.0, 3.0, 4.0};
    double time = 0.5;
    double left_point = 0.0;
    double right_point = 5.0;

    save_step(result, index, y, time, left_point, right_point);

    TEST_ASSERT_EQUAL_DOUBLE(left_point, result->solution_y[index][0]);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(y, &result->solution_y[index][1], computation_grid->N);
    TEST_ASSERT_EQUAL_DOUBLE(right_point, result->solution_y[index][result->grid_size - 1]);
    TEST_ASSERT_EQUAL_DOUBLE(time, result->solution_time[index]);

    // Don't forget to clean up
    destroy_return_data(result);
    destroy_grid(computation_grid);
}

TEST_GROUP_RUNNER(ReturnDataTests)
{
    RUN_TEST_CASE(ReturnDataTests, create_and_destroy_grid);
    RUN_TEST_CASE(ReturnDataTests, save_step_basic);
}

void ReturnDataTests(void)
{
    RUN_TEST_GROUP(ReturnDataTests);
}
