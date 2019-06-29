# fast-math and SIMD instruction settings below has been copied and modified from 
# th GLM library CMakeLists.txt (MIT license)
#
# https://github.com/g-truc/glm/blob/master/CMakeLists.txt

option(ASMC_ENABLE_FAST_MATH "Enable fast math optimizations" OFF)
if(ASMC_ENABLE_FAST_MATH)
  message(STATUS "ASMC: Build with fast math optimizations")

	if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
		add_compile_options(-ffast-math)

	elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		add_compile_options(/fp:fast)
	endif()
else()
	if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		add_compile_options(/fp:precise)
	endif()
endif()

option(ASMC_ENABLE_SIMD_SSE2 "Enable SSE2 optimizations" OFF)
option(ASMC_ENABLE_SIMD_SSE3 "Enable SSE3 optimizations" OFF)
option(ASMC_ENABLE_SIMD_SSSE3 "Enable SSSE3 optimizations" OFF)
option(ASMC_ENABLE_SIMD_SSE4_1 "Enable SSE 4.1 optimizations" OFF)
option(ASMC_ENABLE_SIMD_SSE4_2 "Enable SSE 4.2 optimizations" OFF)
option(ASMC_ENABLE_SIMD_AVX "Enable AVX optimizations" OFF)
option(ASMC_FORCE_PURE "Force 'pure' instructions" OFF)

if(ASMC_FORCE_PURE)
	add_definitions(-DNO_SSE)

	if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
		add_compile_options(-mfpmath=387)
	endif()
	message(STATUS "ASMC: No SIMD instruction set")

elseif(ASMC_ENABLE_SIMD_AVX)
  add_definitions(-DAVX)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-mavx)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxAVX)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		add_compile_options(/arch:AVX)
	endif()
	message(STATUS "ASMC: AVX instruction set")

elseif(ASMC_ENABLE_SIMD_SSE4_2)
	add_definitions(-DGLM_FORCE_INTRINSICS)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-msse4.2)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxSSE4.2)
	elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
		add_compile_options(/arch:SSE2) # VC doesn't support SSE4.2
	endif()
	message(STATUS "ASMC: SSE4.2 instruction set")

elseif(ASMC_ENABLE_SIMD_SSE4_1)
	add_definitions(-DGLM_FORCE_INTRINSICS)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-msse4.1)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxSSE4.1)
	elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
		add_compile_options(/arch:SSE2) # VC doesn't support SSE4.1
	endif()
	message(STATUS "ASMC: SSE4.1 instruction set")

elseif(ASMC_ENABLE_SIMD_SSSE3)
	add_definitions(-DGLM_FORCE_INTRINSICS)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-mssse3)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxSSSE3)
	elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
		add_compile_options(/arch:SSE2) # VC doesn't support SSSE3
	endif()
	message(STATUS "ASMC: SSSE3 instruction set")

elseif(ASMC_ENABLE_SIMD_SSE3)
	add_definitions(-DGLM_FORCE_INTRINSICS)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-msse3)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxSSE3)
	elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
		add_compile_options(/arch:SSE2) # VC doesn't support SSE3
	endif()
	message(STATUS "ASMC: SSE3 instruction set")

elseif(ASMC_ENABLE_SIMD_SSE2)
	add_definitions(-DGLM_FORCE_INTRINSICS)

	if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		add_compile_options(-msse2)
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		add_compile_options(/QxSSE2)
	elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
		add_compile_options(/arch:SSE2)
	endif()
	message(STATUS "ASMC: SSE2 instruction set")
endif()

