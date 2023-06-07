set(LIBRARY_SOURCES 
  grid.c
  return_data.c
  computation_data.c
  computation_library.c
  helpers.c
)

add_library(sim STATIC ${LIBRARY_SOURCES})
target_link_libraries(sim SUNDIALS::cvode m)
target_include_directories(sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB MODELS ${CMAKE_CURRENT_SOURCE_DIR}/models/*)
foreach(FILE ${MODELS})
  get_filename_component(FILENAME ${FILE} NAME)
  get_filename_component(LIBNAME ${FILE} NAME_WE)
  message("File: ${FILENAME}")
  message("Libname: lib${LIBNAME}")
  add_library(${LIBNAME} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/models/${FILENAME} ${LIBRARY_SOURCES})
  target_link_libraries(${LIBNAME} SUNDIALS::cvode m)
  target_include_directories(${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if (${LOG_LINE})
      target_compile_definitions(${LIBNAME} PUBLIC LOG_LINE)
  endif()
endforeach()
