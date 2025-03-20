set(MODS_CONFIGURATION  AngularData.cpp
                        ConfigGenerator.cpp
                        ElectronInfo.cpp
                        HamiltonianMatrix.cpp
                        Level.cpp
                        LevelMap.cpp
                        NonRelConfiguration.cpp
                        Projection.cpp
                        RelativisticConfiguration.cpp
                        RelativisticConfigList.cpp
                        CACHE INTERNAL "")
add_library(configuration "${MODS_CONFIGURATION}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(configuration PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(configuration PUBLIC OpenMP::OpenMP_CXX)
endif()

if(USE_KOKKOS)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_KOKKOS)
  find_package(Kokkos 4.2 REQUIRED CONFIG)
  target_link_libraries(configuration PUBLIC Kokkos::kokkos)
  target_include_directories(configuration PRIVATE ${Kokkos_INCLDUE_DIRS})
endif()

target_link_libraries(configuration PUBLIC Eigen3::Eigen)
target_link_libraries(configuration PUBLIC absl::flat_hash_map)
target_include_directories(configuration PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(configuration PUBLIC mbpt)
