# This file is part of ASMC, developed by Pier Francesco Palamara.
#
# ASMC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ASMC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.12)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(ASMC LANGUAGES CXX)

# ensure C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#export compiler flags for code completion engines
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# Set global compiler warnings
if(MSVC)
    add_compile_options(/W3 /WX)
else()
    add_compile_options(-Wall -pedantic -fPIC)
endif()

# Set SIMD compiler options
include(cmake/SIMD.cmake)

set(ASMC_source_dir ${ASMC_SOURCE_DIR}/ASMC_SRC/SRC)
set(ASMC_test_dir ${ASMC_SOURCE_DIR}/ASMC_SRC/TESTS)
set(ASMC_file_dir ${ASMC_SOURCE_DIR}/FILES)
add_definitions(-DASMC_FILE_DIR=\"${ASMC_file_dir}\")

set(ASMC_source
  ${ASMC_source_dir}/HMM.cpp
  ${ASMC_source_dir}/Data.cpp
  ${ASMC_source_dir}/DecodingParams.cpp
  ${ASMC_source_dir}/DecodingQuantities.cpp
  ${ASMC_source_dir}/FileUtils.cpp
  ${ASMC_source_dir}/Individual.cpp
  ${ASMC_source_dir}/MemoryUtils.cpp
  ${ASMC_source_dir}/StringUtils.cpp
  ${ASMC_source_dir}/Timer.cpp
  ${ASMC_source_dir}/ASMC.cpp
  )

set(ASMC_headers
  ${ASMC_source_dir}/HMM.hpp
  ${ASMC_source_dir}/Data.hpp
  ${ASMC_source_dir}/DecodingParams.hpp
  ${ASMC_source_dir}/DecodingQuantities.hpp
  ${ASMC_source_dir}/FileUtils.hpp
  ${ASMC_source_dir}/Individual.hpp
  ${ASMC_source_dir}/MemoryUtils.hpp
  ${ASMC_source_dir}/StringUtils.hpp
  ${ASMC_source_dir}/Timer.hpp
  ${ASMC_source_dir}/ASMC.hpp
  )

set(ASMC_main
  ${ASMC_source_dir}/main.cpp
  )

set(ASMC_unit_tests
  ${ASMC_test_dir}/unit_tests.cpp
  ${ASMC_test_dir}/test_HMM.cpp
  ${ASMC_test_dir}/test_decoding_params.cpp
  )

set(ASMC_regression
  ${ASMC_test_dir}/test_regression.cpp
  )

option(ASMC_PYTHON_BINDINGS "Whether to build the ASMC python bindings" OFF)
if (ASMC_PYTHON_BINDINGS)
    add_subdirectory(pybind11)
    pybind11_add_module(pyASMC ${ASMC_source_dir}/pybind.cpp)
    target_link_libraries(pyASMC PRIVATE ASMC)
endif ()

# The main ASMC library
add_library(ASMC ${ASMC_source} ${ASMC_headers})
target_include_directories(ASMC PUBLIC ${ASMC_source_dir})

# The main ASMC executable
add_executable(ASMC_exe ${ASMC_main})
target_link_libraries(ASMC_exe PRIVATE ASMC)

enable_testing()

add_executable(ASMC_unit_tests ${ASMC_unit_tests})
target_include_directories(ASMC_unit_tests PRIVATE ${ASMC_test_dir})
target_link_libraries(ASMC_unit_tests PRIVATE ASMC)
add_test(test1 ASMC_unit_tests)

add_executable(ASMC_regression ${ASMC_regression})
target_include_directories(ASMC_regression PRIVATE ${ASMC_test_dir})
target_link_libraries(ASMC_regression PRIVATE ASMC)
add_test(regression ASMC_regression)

# boost as a required dependency for ASMC
find_package(Boost REQUIRED COMPONENTS program_options iostreams)
target_link_libraries(ASMC PRIVATE ${Boost_LIBRARIES})
target_include_directories(ASMC PUBLIC ${Boost_INCLUDE_DIR})

# Eigen is required
find_package(Eigen3 REQUIRED)
target_include_directories(ASMC PUBLIC ${EIGEN3_INCLUDE_DIRS})

# zlib is required (at least on Linux) for boost iostreams
find_package(ZLIB REQUIRED)
target_link_libraries(ASMC PRIVATE ZLIB::ZLIB)
target_include_directories(ASMC PUBLIC ${ZLIB_INCLUDE_DIRS})

# openmp is a required dependency
find_package(OpenMP REQUIRED)
add_definitions(-DHAVE_OPENMP)
target_link_libraries(ASMC PRIVATE OpenMP::OpenMP_CXX)


# Optionally run regression tests to ensure core functionality remains unchanged
option(ASMC_REGRESSION_TESTING "Whether to generate regression testing target" OFF)
if (ASMC_REGRESSION_TESTING)

    # Python is required for running regression tests
    find_package(Python3 REQUIRED)

    add_custom_target(
            Regression
            COMMAND ${Python3_EXECUTABLE} regression_test.py $<TARGET_FILE:ASMC_exe>
            DEPENDS ASMC_exe
            COMMENT "Running ASMC regression test"
            WORKING_DIRECTORY ${ASMC_test_dir}
    )

endif ()


# Clang tidy as optional static analyzer
option(Template_USE_CLANG_TIDY "Use clang tidy for static analysis" OFF)
if (Template_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-6.0" "clang-tidy-5.0" "clang-tidy-4.0"
                 DOC "Path to clang-tidy executable")
    if(NOT CLANG_TIDY_EXE)
      message(STATUS "clang-tidy not found.")
    else()
      message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
      set_target_properties(ASMC PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif()
endif()

# LLVM AddressSanitizer to detecting memory errors
# Note that there are many other sanitizers in LLVM to help detect errors, see
# http://travistoptips.blogspot.co.uk/2015/11/sanitize-all-things.html
option(Template_MEMCHECK "Use LLVM AddressSanitizer to detecting memory errors" OFF)
if (Template_MEMCHECK)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Configuring with LLVM AddressSanitizer")
        set(Template_MEMCHECK_FLAGS -fno-optimize-sibling-calls 
                                    -fsanitize=address 
                                    -fsanitize-address-use-after-scope
                                    )
        target_compile_options(ASMC PUBLIC -O1 -g -fno-omit-frame-pointer ${Template_MEMCHECK_FLAGS})
        target_link_libraries(ASMC PUBLIC -g ${Template_MEMCHECK_FLAGS})
    else()
        message(FATAL_ERROR "clang compiler required with Template_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

# Setup coverage testing for GCC or Clang
option(ASMC_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (ASMC_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message(STATUS "Configuring with coverage")
        target_compile_options(ASMC PUBLIC --coverage -O0)
        target_link_libraries(ASMC PUBLIC --coverage)
    else()
        message(FATAL_ERROR "GCC or Clang required with ASMC_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()
