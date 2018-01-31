cmake_minimum_required(VERSION 2.8.12)
project(ppproject_native)

# Detect C++ version support.
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
  if(COMPILER_SUPPORTS_CXX14)
    message("Enabling C++14 support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    set(TARGET_CXX_FLAGS "${TARGET_CXX_FLAGS} -std=c++14")
  elseif(COMPILER_SUPPORTS_CXX14)
    message(WARNING "No C++14 support, enabling C++11 support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(TARGET_CXX_FLAGS "${TARGET_CXX_FLAGS} -std=c++11")
  else()
    message(WARNING "No C++11 support. This may cause issues.")
  endif()
endif()

# Require pybind11.
#find_package(pybind11)
add_subdirectory("../../dep/pybind11" "pybind11")

# Set native folder as an include
include_directories(${CMAKE_SOURCE_DIR})

# Shared library
add_subdirectory(shared)

# Modules
pybind11_add_module(samplebuffer samplebuffer.cpp)
target_link_libraries(samplebuffer PRIVATE shared)

