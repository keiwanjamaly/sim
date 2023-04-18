add_executable(test_runner all_tests.c test_grid.c)

# Link libraries
target_link_libraries(test_runner PRIVATE grid_code Unity::Unity)

# Add the UNITY_EXTENSION_FIXTURE definition
target_compile_definitions(test_runner PRIVATE UNITY_EXTENSION_FIXTURE=1)

# Add a custom command to run the tests after building
add_custom_command(TARGET test_runner POST_BUILD COMMAND test_runner)
