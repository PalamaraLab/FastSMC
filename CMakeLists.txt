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
include(cmake/AutodetectVcpkgToolchainFile.cmake)

project(asmc-fastsmc LANGUAGES CXX)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    message(STATUS "  You can specify with the flag -DCMAKE_BUILD_TYPE=<Debug|Release|MinSizeRel|RelWithDebInfo>")
    message(STATUS "  'Release' will build optimised binaries, but 'Debug' may be better while developing FastSMC.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()

message(STATUS "CMake build type is set to ${CMAKE_BUILD_TYPE}")

# ensure C++17
set(CMAKE_CXX_STANDARD 17)
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

set(ASMC_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/ASMC_SRC/SRC)
set(ASMC_test_dir ${CMAKE_CURRENT_SOURCE_DIR}/ASMC_SRC/TESTS)
set(ASMC_file_dir ${CMAKE_CURRENT_SOURCE_DIR}/FILES)
add_definitions(-DASMC_TEST_DIR=\"${ASMC_test_dir}\")
add_definitions(-DASMC_FILE_DIR=\"${ASMC_file_dir}\")

set(ASMC_source
  ${ASMC_source_dir}/HMM.cpp
  ${ASMC_source_dir}/HmmUtils.cpp
  ${ASMC_source_dir}/Data.cpp
  ${ASMC_source_dir}/DecodingParams.cpp
  ${ASMC_source_dir}/DecodingQuantities.cpp
  ${ASMC_source_dir}/FastSMC.cpp
  ${ASMC_source_dir}/FileUtils.cpp
  ${ASMC_source_dir}/Individual.cpp
  ${ASMC_source_dir}/MemoryUtils.cpp
  ${ASMC_source_dir}/StringUtils.cpp
  ${ASMC_source_dir}/Timer.cpp
  ${ASMC_source_dir}/ASMC.cpp
  ${ASMC_source_dir}/HASHING/Utils.cpp
  )

set(ASMC_headers
  ${ASMC_source_dir}/AvxDefinitions.hpp
  ${ASMC_source_dir}/BinaryDataReader.hpp
  ${ASMC_source_dir}/HMM.hpp
  ${ASMC_source_dir}/HmmUtils.hpp
  ${ASMC_source_dir}/Data.hpp
  ${ASMC_source_dir}/DecodePairsReturnStruct.hpp
  ${ASMC_source_dir}/DecodingParams.hpp
  ${ASMC_source_dir}/DecodingQuantities.hpp
  ${ASMC_source_dir}/FastSMC.hpp
  ${ASMC_source_dir}/FileUtils.hpp
  ${ASMC_source_dir}/Individual.hpp
  ${ASMC_source_dir}/MemoryUtils.hpp
  ${ASMC_source_dir}/StringUtils.hpp
  ${ASMC_source_dir}/Timer.hpp
  ${ASMC_source_dir}/ASMC.hpp
  ${ASMC_source_dir}/HASHING/ExtendHash.hpp
  ${ASMC_source_dir}/HASHING/Individuals.hpp
  ${ASMC_source_dir}/HASHING/Match.hpp
  ${ASMC_source_dir}/HASHING/SeedHash.hpp
  ${ASMC_source_dir}/HASHING/Utils.hpp
  )

set(ASMC_main
  ${ASMC_source_dir}/main.cpp
  )

set(FastSMC_main
  ${ASMC_source_dir}/main_fastsmc.cpp
  )

set(convertBinary_main
  ${ASMC_source_dir}/main_convertBinary.cpp
  )

set(ASMC_unit_tests
  ${ASMC_test_dir}/unit_tests.cpp
  ${ASMC_test_dir}/test_ASMC.cpp
  ${ASMC_test_dir}/test_binary_data_reader.cpp
  ${ASMC_test_dir}/test_hashing.cpp
  ${ASMC_test_dir}/test_HMM.cpp
  ${ASMC_test_dir}/test_hmm_utils.cpp
  ${ASMC_test_dir}/test_decoding_params.cpp
  ${ASMC_test_dir}/test_decoding_quantities.cpp
  )

set(ASMC_regression
  ${ASMC_test_dir}/test_regression.cpp
  ${ASMC_test_dir}/test_fastsmc_regression.cpp
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

# Link against std filesystem on GCC < 8.4
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.4)
        target_link_libraries(ASMC PRIVATE stdc++fs)
    endif ()
endif ()

# The main ASMC executable
add_executable(ASMC_exe ${ASMC_main})
target_link_libraries(ASMC_exe PRIVATE ASMC)

# The main FastSMC executable
add_executable(FastSMC_exe ${FastSMC_main})
target_link_libraries(FastSMC_exe PRIVATE ASMC)

# The main convertBinary executable
add_executable(convertBinary_exe ${convertBinary_main})
target_link_libraries(convertBinary_exe PRIVATE ASMC)

enable_testing()

add_executable(ASMC_unit_tests ${ASMC_unit_tests})
target_include_directories(ASMC_unit_tests PRIVATE ${ASMC_test_dir})
target_link_libraries(ASMC_unit_tests PRIVATE ASMC)
add_test(Asmc_unit_tests ASMC_unit_tests)

add_executable(ASMC_regression ${ASMC_regression})
target_include_directories(ASMC_regression PRIVATE ${ASMC_test_dir})
target_link_libraries(ASMC_regression PRIVATE ASMC)
add_test(regression ASMC_regression)

# boost as a required dependency for ASMC
find_package(Boost REQUIRED COMPONENTS program_options iostreams)
target_link_libraries(ASMC PRIVATE ${Boost_LIBRARIES})
target_include_directories(ASMC PUBLIC ${Boost_INCLUDE_DIR})

# Eigen is required
find_package(Eigen3 CONFIG REQUIRED)
target_include_directories(ASMC PUBLIC ${EIGEN3_INCLUDE_DIRS})

# zlib is required (at least on Linux) for boost iostreams
find_package(ZLIB REQUIRED)
target_link_libraries(ASMC PRIVATE ZLIB::ZLIB)
target_include_directories(ASMC PUBLIC ${ZLIB_INCLUDE_DIRS})

# openmp is a required dependency
find_package(OpenMP REQUIRED)
add_definitions(-DHAVE_OPENMP)
target_link_libraries(ASMC PRIVATE OpenMP::OpenMP_CXX)

# fmt is required for efficient string formatting
find_package(fmt CONFIG REQUIRED)
target_link_libraries(ASMC PRIVATE fmt::fmt)

# requires asmc-prepare-decoding
#find_package(AsmcDataModule CONFIG REQUIRED)
#target_link_libraries(ASMC PRIVATE AsmcDataModule::data_module_lib)


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
option(ASMC_MEMCHECK "Use LLVM AddressSanitizer to detecting memory errors" OFF)
if (ASMC_MEMCHECK)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Configuring with LLVM AddressSanitizer")
        set(ASMC_MEMCHECK_FLAGS -fno-optimize-sibling-calls
                                    -fsanitize=address 
                                    -fsanitize-address-use-after-scope
                                    )
        target_compile_options(ASMC PUBLIC -O1 -g -fno-omit-frame-pointer ${ASMC_MEMCHECK_FLAGS})
        target_link_libraries(ASMC PUBLIC -g ${ASMC_MEMCHECK_FLAGS})
    else()
        message(FATAL_ERROR "clang compiler required with ASMC_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
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
