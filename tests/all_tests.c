#include "unity_fixture.h"

static void RunAllTests(void)
{
    // Add calls to RUN_TEST_GROUP() for each test group here
    RUN_TEST_GROUP(GridTests);
}

int main(int argc, char *argv[])
{
    return UnityMain(argc, argv, RunAllTests);
}