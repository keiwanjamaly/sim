set(LIBRARY_SOURCES 
  grid.c
  return_data.c
  computation_data.c
)

add_library(libsim STATIC ${LIBRARY_SOURCES})
target_include_directories(libsim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)