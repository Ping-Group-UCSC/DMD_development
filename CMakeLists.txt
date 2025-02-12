cmake_minimum_required(VERSION 3.10)
project(denmat_dynm_v4.5.7 CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optionally set build type if not specified externally
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -traceback")
endif()

# --- Compiler Flags ---
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# --- GSL Setup ---
set(GSL_INCLUDE_DIR "/usr/include/gsl")
if(EXISTS "${GSL_INCLUDE_DIR}")
    message(STATUS "Found GSL include directory: ${GSL_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "GSL include directory not found: ${GSL_INCLUDE_DIR}")
endif()

# --- Option to Use MKL ---
option(USE_MKL "Use MKL for linear algebra operations" ON)

if(USE_MKL AND DEFINED ENV{MKLROOT})
    # --- MKL Setup ---
    set(MKLROOT $ENV{MKLROOT})
    message(STATUS "Using MKL from MKLROOT: ${MKLROOT}")

    # Locate the MKL include directory.
    if(EXISTS "${MKLROOT}/latest/include")
        set(MKL_INCLUDE_DIR "${MKLROOT}/latest/include")
    elseif(EXISTS "${MKLROOT}/include")
        set(MKL_INCLUDE_DIR "${MKLROOT}/include")
    else()
        message(FATAL_ERROR "MKL include directory not found in ${MKLROOT}/latest/include or ${MKLROOT}/include")
    endif()

    # Locate the MKL library directory.
    if(EXISTS "${MKLROOT}/latest/lib/intel64")
        set(MKL_LIB_DIR "${MKLROOT}/latest/lib/intel64")
    elseif(EXISTS "${MKLROOT}/lib/intel64")
        set(MKL_LIB_DIR "${MKLROOT}/lib/intel64")
    else()
        message(FATAL_ERROR "MKL library directory not found in ${MKLROOT}/latest/lib/intel64 or ${MKLROOT}/lib/intel64")
    endif()
    message(STATUS "MKL include directory: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL library directory: ${MKL_LIB_DIR}")
    link_directories(${MKL_LIB_DIR})

    # Define core MKL linker flags (for the LP64 interface using sequential threading).
    set(MKL_LIBS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")

    # --- MKL Scalapack and BLACS Setup (for MPI) ---
    # Adjust these if you are using a different MPI implementation or integer interface.
    set(MKL_SCALAPACK_LIB "mkl_scalapack_lp64")
    set(MKL_BLACS_LIB "mkl_blacs_openmpi_lp64")
else()
    message(STATUS "MKL not used or MKLROOT not defined. Falling back to system Scalapack.")
    set(USE_MKL OFF)
    # --- System SCALAPACK Setup ---
    # These paths are typical on Linux; adjust them for macOS (e.g., if installed via Homebrew).
    set(SCALAPACK_INCLUDE_DIR "/usr/include")
    set(SCALAPACK_LIB_DIR "/usr/lib/x86_64-linux-gnu")
    message(STATUS "Assuming system SCALAPACK headers are in: ${SCALAPACK_INCLUDE_DIR}")
    message(STATUS "Assuming system SCALAPACK libraries are in: ${SCALAPACK_LIB_DIR}")
    # Uncomment the following line if needed:
    # link_directories(${SCALAPACK_LIB_DIR})
endif()

# --- Project Source Files ---
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/src_v4.5.7/*.cpp")
if(NOT SOURCES)
    message(FATAL_ERROR "No source files found in ${CMAKE_SOURCE_DIR}/src_v4.5.7")
endif()
message(STATUS "Source files found:")
foreach(src ${SOURCES})
    message(STATUS "  ${src}")
endforeach()

# --- Include Directories for the Project ---
set(PROJECT_INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}/src_v4.5.7"
    "${CMAKE_SOURCE_DIR}/src_v4.5.7/comm"
    "${CMAKE_SOURCE_DIR}/src_v4.5.7/material_model"
    "${CMAKE_SOURCE_DIR}/src_v4.5.7/scattering_model"
    "${CMAKE_SOURCE_DIR}/src_v4.5.7/material_model/GaAs"
    "${CMAKE_SOURCE_DIR}/src_v4.5.7/material_model/MoS2"
)
foreach(dir ${PROJECT_INCLUDE_DIRS})
    if(EXISTS "${dir}")
        message(STATUS "Found project include directory: ${dir}")
    else()
        message(WARNING "Project include directory not found: ${dir}")
    endif()
endforeach()

# --- Create the Executable ---
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Target-Specific Include Directories ---
target_include_directories(${PROJECT_NAME} PRIVATE
    ${PROJECT_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIR}
)

if(USE_MKL)
    target_include_directories(${PROJECT_NAME} PRIVATE ${MKL_INCLUDE_DIR})
else()
    target_include_directories(${PROJECT_NAME} PRIVATE ${SCALAPACK_INCLUDE_DIR})
endif()

# --- MPI Setup ---
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    message(STATUS "Found MPI C++ compiler: ${MPI_CXX_COMPILER}")
    message(STATUS "MPI include path: ${MPI_CXX_INCLUDE_PATH}")
else()
    message(FATAL_ERROR "MPI C++ compiler not found!")
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH})

# --- Link Libraries ---
if(USE_MKL)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        gsl gslcblas
        ${MKL_LIBS}
        ${MKL_SCALAPACK_LIB}
        ${MKL_BLACS_LIB}
        MPI::MPI_CXX
    )
else()
    target_link_libraries(${PROJECT_NAME} PRIVATE
        gsl gslcblas
        scalapack
        MPI::MPI_CXX
    )
endif()

# --- Set Output Directories ---
set_target_properties(${PROJECT_NAME} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
)
message(STATUS "Executable will be output to: ${CMAKE_SOURCE_DIR}/bin")
