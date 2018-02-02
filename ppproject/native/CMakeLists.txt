cmake_minimum_required(VERSION 2.8.12)
project(ppproject_native)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

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

# Enable PIC globally so that we can link static libs against our python modules.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Shared library
add_subdirectory(shared)

# samplebuffer module
pybind11_add_module(samplebuffer MODULE samplebuffer.cpp)
target_link_libraries(samplebuffer PRIVATE shared)

# filereader module
find_package(LibSndFile REQUIRED)
pybind11_add_module(filereader MODULE filereader.cpp)
target_include_directories(filereader PRIVATE ${LIBSNDFILE_INCLUDE_DIR})
target_link_libraries(filereader PRIVATE shared ${LIBSNDFILE_LIBRARY})

# spectrogram module
find_package(Cairo REQUIRED)
find_package(FFTW REQUIRED)
pybind11_add_module(spectrogram MODULE spectrogram.cpp)
target_include_directories(spectrogram PRIVATE ${CAIRO_INCLUDE_DIRS} ${FFTW_INCLUDES})
target_link_libraries(spectrogram PRIVATE shared ${CAIRO_LIBRARIES} ${FFTW_LIBRARIES})

