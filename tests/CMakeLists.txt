add_executable(test_runner all_tests.c test_grid.c)

# Add Unity header directories
target_include_directories(test_runner PRIVATE
    /usr/local/include/unity
    /usr/local/include/unity/extras/fixture/src
    /usr/local/include/unity/extras/memory/src
)

# Add the UNITY_EXTENSION_FIXTURE and UNITY_DOUBLE_PRECISION definitions
target_compile_definitions(test_runner PRIVATE UNITY_EXTENSION_FIXTURE=1)

# Link libraries
target_link_libraries(test_runner PRIVATE grid_code /usr/local/lib/libunity.a)

# Add a custom command to run the tests after building
add_custom_command(TARGET test_runner POST_BUILD COMMAND test_runner)
