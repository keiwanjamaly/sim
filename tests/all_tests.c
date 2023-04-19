#include "unity_fixture.h"

static void RunAllTests(void)
{
    // Add calls to RUN_TEST_GROUP() for each test group here
    RUN_TEST_GROUP(GridTests);
    RUN_TEST_GROUP(ReturnDataTests);
}

int main(int argc, const char *argv[])
{
    return UnityMain(argc, argv, RunAllTests);
}