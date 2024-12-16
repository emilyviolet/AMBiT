set(MODS_SPECIFICATION Specification.cpp
                    CACHE INTERNAL "")

add_library(specification "${MODS_SPECIFICATION}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(specification PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(specification PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(specification PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(specification PUBLIC OpenMP::OpenMP_CXX)
endif()

