if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  
  # Boilerplate for setting up googletest
  # Download gtest from GitHub
  include(FetchContent)
    FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip)
  FetchContent_MakeAvailable(googletest)

  # Actual ambit testing libraries
  set(MODS_TESTING AngularData.test.cpp
                   BSplineBasis.test.cpp
                   BasisGenerator.test.cpp
                   Breit.test.cpp
                   BruecknerDecorator.test.cpp
                   ConfigGenerator.test.cpp
                   ConfigurationParser.test.cpp
                   CoreValenceIntegrals.test.cpp
                   EJOperator.test.cpp
                   HFOperator.test.cpp
                   HamiltonianMatrix.test.cpp
                   HartreeFocker.test.cpp
                   Hyperfine.test.cpp
                   LocalPotentialDecorator.test.cpp
                   ManyBodyOperator.test.cpp
                   MassShiftDecorator.test.cpp
                   MathConstant.test.cpp
                   MultirunOptions.test.cpp
                   RadiativePotential.test.cpp
                   ambit.test.cpp
                   CACHE INTERNAL "")

  # Build the test exe and register it with gtest
  add_executable(ambit_test "${MODS_TESTING}")

  # Internal AMBiT link dependencies
  target_link_libraries(ambit_test PUBLIC atom)
  target_link_libraries(ambit_test PUBLIC basis)
  target_link_libraries(ambit_test PUBLIC configuration)
  target_link_libraries(ambit_test PUBLIC externalfield)
  target_link_libraries(ambit_test PUBLIC hartreefock)
  target_link_libraries(ambit_test PUBLIC mbpt)
  target_link_libraries(ambit_test PUBLIC universal)

  # External dependencies
  find_package(Eigen3 REQUIRED NO_MODULE)
  target_link_libraries(ambit_test PRIVATE Eigen3::Eigen)
  target_link_libraries(ambit_test PUBLIC GSL::gsl)
  find_package(GSL REQUIRED)

  if(USE_MKL)
    find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})  
    target_link_libraries(ambit_test PUBLIC MKL::MKL)
    # Also tell Eigen to use MKL as its backend
    target_compile_definitions(ambit_test PRIVATE -DEIGEN_USE_MKL_ALL)
  else() # Use vanilla BLAS/LAPACK
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    target_link_libraries(ambit_test PUBLIC LAPACK::LAPACK)
    target_link_libraries(ambit_test PUBLIC BLAS::BLAS)
  endif()

  if(USE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(ambit_test PRIVATE MPI::MPI_CXX)
  endif()
  if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(ambit_test PRIVATE OpenMP::OpenMP_CXX)
  endif()

  find_package(Boost REQUIRED)
  target_include_directories(ambit_test PRIVATE ${Boost_INCLUDE_DIRS})

  # gtest dependencies
  target_link_libraries(ambit_test PUBLIC GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(ambit_test)
endif()
