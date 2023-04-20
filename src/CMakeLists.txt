set(LIBRARY_SOURCES 
  grid.c
  return_data.c
  computation_data.c
  computation_library.c
  ${MODEL_FILE}
)

add_library(sim SHARED ${LIBRARY_SOURCES})
target_link_libraries(sim SUNDIALS::cvode m)
target_include_directories(sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)