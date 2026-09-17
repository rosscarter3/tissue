cmake_minimum_required(VERSION 3.20)
project(tissue2 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TISSUE_NATIVE "Optimize for the build machine's CPU (SIMD)" ON)
option(TISSUE_OPENMP "Use OpenMP if available (falls back to built-in thread pool)" OFF)

add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(TISSUE_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-mcpu=native" HAS_MCPU_NATIVE)
  check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
  if(HAS_MCPU_NATIVE)
    add_compile_options(-mcpu=native)
  elseif(HAS_MARCH_NATIVE)
    add_compile_options(-march=native)
  endif()
endif()

# OBJECT library: keeps the self-registering reaction/solver translation
# units alive (a static archive would let the linker drop them).
file(GLOB_RECURSE TISSUE2_SOURCES CONFIGURE_DEPENDS src/tissue/*.cpp)
add_library(tissue2 OBJECT ${TISSUE2_SOURCES})
target_include_directories(tissue2 PUBLIC src)

find_package(Threads REQUIRED)
target_link_libraries(tissue2 PUBLIC Threads::Threads)

if(TISSUE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(tissue2 PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(tissue2 PUBLIC TISSUE_HAVE_OPENMP=1)
  endif()
endif()

add_executable(simulator src/apps/simulator.cpp)
target_link_libraries(simulator PRIVATE tissue2)

enable_testing()
add_subdirectory(tests)
