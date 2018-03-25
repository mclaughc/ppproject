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
  elseif(COMPILER_SUPPORTS_CXX14)
    message(WARNING "No C++14 support, enabling C++11 support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
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

# External dependencies
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/msvc ${CMAKE_CURRENT_SOURCE_DIR}/build/dep/msvc)
else()
  find_package(LibSndFile REQUIRED)
  find_package(FFTW REQUIRED)

  # Find and enable OpenMP
  # This is in here because MSVC doesn't support max reduction in OpenMP.
  find_package(OpenMP)
  if (OPENMP_FOUND)
    message(STATUS "OpenMP found, enabling.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message(WARNING "OpenMP not found, performance may be slower as a result.")
  endif()
endif()

# Shared library
add_subdirectory(shared)

# samplebuffer module
pybind11_add_module(samplebuffer MODULE samplebuffer.cpp)
target_link_libraries(samplebuffer PRIVATE shared)

# audiofile module
pybind11_add_module(audiofile MODULE audiofile.cpp)
target_include_directories(audiofile PRIVATE ${LIBSNDFILE_INCLUDE_DIR})
target_link_libraries(audiofile PRIVATE shared ${LIBSNDFILE_LIBRARY})

# spectrogram module
pybind11_add_module(spectrogram MODULE spectrogram.cpp)
target_include_directories(spectrogram PRIVATE ${FFTW_INCLUDES})
target_link_libraries(spectrogram PRIVATE shared ${FFTW_LIBRARIES})

# lowpassfilter module
pybind11_add_module(lowpassfilter MODULE lowpassfilter.cpp)
target_link_libraries(lowpassfilter PRIVATE shared)

# copy targets to native directory, so it can be imported from python.
install(TARGETS samplebuffer audiofile spectrogram lowpassfilter DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

# copy shared .dlls to native directory, for win32
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    FILE(GLOB DEP_BINARIES "${CMAKE_CURRENT_SOURCE_DIR}/../../dep/msvc/bin64/*.dll")
  else()
    FILE(GLOB DEP_BINARIES "${CMAKE_CURRENT_SOURCE_DIR}/../../dep/msvc/bin/*.dll")
  endif()
  install(FILES ${DEP_BINARIES} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
endif()
