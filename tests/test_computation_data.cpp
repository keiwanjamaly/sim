#include "unity.h"
#include "unity_fixture.h"
#include "computation_data.h"

struct physics_data
{
    int i;
};

TEST_GROUP(ComputationDataTests);

TEST_SETUP(ComputationDataTests)
{
}

TEST_TEAR_DOWN(ComputationDataTests)
{
}

TEST(ComputationDataTests, create_and_destroy_grid)
{
    double Lambda = 1.0;
    double kir = 0.5;
    Grid computation_grid;
    PhysicsData data;
    double tolerances = 0.01;

    ComputationData *result = initialize_computation_data(Lambda, kir, &computation_grid, &data, tolerances);

    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_DOUBLE(Lambda, result->Lambda);
    TEST_ASSERT_EQUAL_DOUBLE(kir, result->kir);
    TEST_ASSERT_EQUAL_DOUBLE(-log(kir / Lambda), result->tir);
    TEST_ASSERT_EQUAL_DOUBLE(tolerances, result->tolerances);
    TEST_ASSERT_EQUAL_PTR(&computation_grid, result->computation_grid);
    TEST_ASSERT_EQUAL_PTR(&data, result->data);

    destroy_computation_data(result);
}

TEST_GROUP_RUNNER(ComputationDataTests)
{
    RUN_TEST_CASE(ComputationDataTests, create_and_destroy_grid);
}

void ComputationDataTests(void)
{
    RUN_TEST_GROUP(ComputationDataTests);
}
