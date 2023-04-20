# Add your test files here
set(TEST_SOURCES
  all_tests.c
  test_grid.c
  test_return_data.c
  test_computation_data.c
)

# Create the test executable
add_executable(run_tests ${TEST_SOURCES})

# target_compile_definitions(run_tests PUBLIC UNITY_INCLUDE_DOUBLE)

# Link Unity to the test executable
target_link_libraries(run_tests unity::framework sim)

# target_include_directories(run_tests unity)
# target_include_directories(run_tests PUBLIC ${UNITY_INCLUDE_DIRS})
