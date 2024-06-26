set(MODS_ATOM Atom.cpp
              Atom_Open.cpp
              Atom_Autoionization.cpp
              OutStreams.cpp
              CACHE INTERNAL "")

add_library(atom "${MODS_ATOM}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(atom PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(atom PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(atom PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(atom PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(atom PUBLIC externalfield)
