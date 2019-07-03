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


cmake_minimum_required(VERSION 3.6)
project(ASMC)

# ensure C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#export compiler flags for code completion engines
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# Set global compiler warnings
if(MSVC)
    add_compile_options(/W3 /WX)
else()
    add_compile_options(-Wall -pedantic)
endif()

# Set SIMD compiler options
include(cmake/SIMD.cmake)

set(ASMC_source_dir ASMC_SRC/SRC)
set(ASMC_test_dir ASMC_SRC/TESTS)
set(ASMC_source
  ${ASMC_source_dir}/Data.cpp
  ${ASMC_source_dir}/DecodingParams.cpp
  ${ASMC_source_dir}/DecodingQuantities.cpp
  ${ASMC_source_dir}/FileUtils.cpp
  ${ASMC_source_dir}/Individual.cpp
  ${ASMC_source_dir}/MemoryUtils.cpp
  ${ASMC_source_dir}/StringUtils.cpp
  ${ASMC_source_dir}/Timer.cpp
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
  )

set(ASMC_main
  ${ASMC_source_dir}/main.cpp
  )


set(ASMC_test_dir ASMC_SRC/TESTS)
set(ASMC_tests
  ${ASMC_test_dir}/Placeholder.cpp
  )


add_library(ASMC ${ASMC_source} ${ASMC_headers})
target_include_directories(ASMC PUBLIC ${ASMC_source_dir})
add_executable(ASMC_exe ${ASMC_main})
target_link_libraries(ASMC_exe PRIVATE ASMC)

enable_testing()

add_executable(ASMC_tests ${ASMC_tests})
target_include_directories(ASMC_tests PRIVATE ${ASMC_test_dir})
target_link_libraries(ASMC_tests PRIVATE ASMC)
add_test(test1 ASMC_tests)

# boost as a required dependency for ASMC
find_package(Boost REQUIRED COMPONENTS program_options iostreams)
target_link_libraries(ASMC PRIVATE "${Boost_LIBRARIES}")
target_include_directories(ASMC PRIVATE "${Boost_INCLUDE_DIR}")

# openmp is a required dependency
find_package(OpenMP REQUIRED)
add_definitions(-DHAVE_OPENMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

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
option(Template_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (Template_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message(STATUS "Configuring with coverage")
        target_compile_options(ASMC PUBLIC --coverage -O0)
        target_link_libraries(ASMC PUBLIC --coverage)
    else()
        message(FATAL_ERROR "GCC or Clang required with Template_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()