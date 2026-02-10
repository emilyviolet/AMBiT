# First, generate the gitinfo header file and dummy target
execute_process(COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/getGitInfo.sh" 
                OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/gitInfo.h")
add_custom_target(generate_gitinfo ALL DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/gitInfo.h")

# Main ambit exe
add_executable(ambit ambit.cpp ambit.h)
# Need to explicitly add the dependency on gitInfo.h to ensure it gets built before we do anything
# else
add_dependencies(ambit generate_gitinfo)

# External libraries
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(absl CONFIG REQUIRED)
target_link_libraries(ambit PRIVATE absl::strings)
target_link_libraries(ambit PRIVATE absl::str_format)

# Use GPU-accelerated libraries?
set(USE_GPU OFF CACHE BOOL "Use GPU-accelerated subroutines")
set(GPU_BACKEND CACHE STRING "Backend API to use for GPU-offloading. Can be CUDA or HIP")
set(GPU_BACKEND_VALUES "HIP" "CUDA")
set_property(CACHE GPU_BACKEND PROPERTY STRINGS ${GPU_BACKEND_VALUES})
string(TOUPPER ${GPU_BACKEND_VALUES} GPU_BACKEND_VALUES)

# LAPACK and BLAS
# First, check to see if we need to use MKL or OpenBLAS. Can't use both at once, so error out if
# both are set.
set(USE_MKL OFF CACHE BOOL "Use Intel MKL for libear algebra subroutines")
set(USE_OPENBLAS OFF CACHE BOOL "Use OpenBLAS for linear algebra subroutines")
set(USE_SCALAPACK OFF CACHE BOOL "Use ScaLAPACK when calculating the CI matrix eigenvalues")

# Verification to make sure there are no conflicting linear algebra libraries
if(USE_MKL AND USE_OPENBLAS)
    message(FATAL_ERROR "Cannot use both MKL and OpenBLAS for linear algebra. Only enable one at a
    time.")
endif()

if(USE_SCALAPACK AND NOT USE_MPI)
    message(FATAL_ERROR "Must enable MPI to use ScaLAPACK.")
endif()

if(USE_MKL)
    # MKL 2021.3 and later has a nice CMake support, so make sure we use that
    # Documentation can be found here:
    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/cmake-config-for-onemkl.html
    # and
    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/using-the-cmake-config-file.html
    # Get the path from the shell or CMake command-line options (e.g. set -DMKL_ROOT=/path/to/MKL)
    find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
    target_link_libraries(ambit PRIVATE MKL::MKL)
    # Also tell Eigen to use MKL as its backend
    target_compile_definitions(ambit PRIVATE -DEIGEN_USE_MKL_ALL)

    # Use the MKL ScaLAPACK backend if requested
    if(USE_SCALAPACK)
        target_compile_definitions(ambit PRIVATE -DAMBIT_USE_SCALAPACK)
        set(ENABLE_SCALAPACK true)
        set(ENABLE_BLACS true)
        target_link_libraries(ambit PRIVATE MKL::MKL_SCALAPACK)
    endif()

else() # Use vanilla BLAS/LAPACK

    if(USE_OPENBLAS)
        # Need to tell Eigen to use external BLAS lib via -DEIGEN_USE_BLAS, then manually pass the
        # openblas linker flags through CMake's "BLA_VENDOR" config variable
        # See this page in the Eigen documentation for an explanation of this process: https://eigen.tuxfamily.org/dox/TopicUsingBlasLapack.html
        set(BLA_VENDOR OpenBLAS)
        target_compile_definitions(ambit PRIVATE -DEIGEN_USE_BLAS)
    endif()
    find_package(LAPACK REQUIRED)

    target_link_libraries(ambit PRIVATE LAPACK::LAPACK)
    find_package(BLAS REQUIRED)
    target_link_libraries(ambit PRIVATE BLAS::BLAS)

    if(USE_SCALAPACK)
        target_compile_definitions(ambit PRIVATE -DAMBIT_USE_SCALAPACK)
        find_package(scalapack REQUIRED)
        target_link_libraries(ambit PRIVATE scalapack)
    endif()

endif()

# GPU options
if(USE_GPU)
    # Choose linking based on whether we're using HIP or CUDA for linear algebra on the
    # GPU
    if(GPU_BACKEND STREQUAL "HIP")
        # AMD libraries
        find_package(rocblas    REQUIRED)
        find_package(rocsolver  REQUIRED)
        target_link_libraries(ambit PRIVATE rocblas)
        target_link_libraries(ambit PRIVATE rocsolver)

        # Set the internal AMBiT preprocessor definition to keep track of which GPU
        # backend we're using (has to be a boolean because of the preprocessor)
        target_compile_definitions(ambit PRIVATE -DAMBIT_GPU_HIP)

    elseif(GPU_BACKEND STREQUAL "CUDA")
        # NVIDIA libraries
        find_package(cublas    REQUIRED)
        find_package(cusolver  REQUIRED)
        target_link_libraries(ambit PRIVATE cublas)
        target_link_libraries(ambit PRIVATE cusolver)

        # Set the internal AMBiT preprocessor definition to keep track of which GPU
        # backend we're using (has to be a boolean because of the preprocessor)
        target_compile_definitions(ambit PRIVATE -DAMBIT_GPU_CUDA)
        
        message(FATAL_ERROR "CUDA GPU backend not yet implemented")

    else()
        # Bail out if passed an unsupported GPU backend
        message(FATAL_ERROR "Unsupported GPU backend ${GPU_BACKEND}")
    endif()

    # Internal AMBIT preprocessor definition to activate GPU compilation paths
    target_compile_definitions(ambit PRIVATE -DAMBIT_USE_GPU)

endif()

# MPI options

# Important note: our user-facing "CACHE" variables (i.e. we can set them with cmake -DUSE_MPI=yes)
# have to be different from the internal AMBiT preprocessor definitions. AMBiT just checks if
# AMBIT_USE_(OPENMP|MPI) are defined, regardless of its actual value. If we expose these variables
# in the CMake cache then they'll get set to FALSE, but the preprocessor will see that they exist
# and thus will always generate the OpenMP/MPI code, even though that's not what the user has
# requested. Consequently, we need "PRIVATE" compile definitions that shadow the "CACHE" ones in
# order to make this work properly.
set(USE_MPI OFF CACHE BOOL "Enable MPI support in AMBiT")
if(USE_MPI)
  target_compile_definitions(ambit PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(ambit PUBLIC MPI::MPI_CXX)
endif()

# OpenMP options. CMake should be able to automatically infer the correct flags for the current
# compiler
set(USE_OPENMP OFF CACHE BOOL "Enable OpenMP support in AMBiT")
if(USE_OPENMP)
  target_compile_definitions(ambit PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(ambit PUBLIC OpenMP::OpenMP_CXX)
endif()

# Now include the directories for the different ambit modules and link against the static libs
add_subdirectory(Atom)
add_subdirectory(Basis)
add_subdirectory(Configuration)
add_subdirectory(ExternalField)
add_subdirectory(HartreeFock)
add_subdirectory(MBPT)
add_subdirectory(Universal)

target_link_libraries(ambit PUBLIC atom basis configuration externalfield hartreefock mbpt universal)

# Install the "ambit" exe into "./bin".
install(TARGETS ambit DESTINATION "bin")
