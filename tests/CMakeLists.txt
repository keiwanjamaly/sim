add_executable(test_runner all_tests.c test_grid.c)
target_link_libraries(test_runner PRIVATE grid_code unity)

# Add a custom command to run the tests after building
add_custom_command(TARGET test_runner POST_BUILD COMMAND test_runner)
