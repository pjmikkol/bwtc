# To compile release version run  'cmake -DCMAKE_BUILD_TYPE=Release'
# For debug version run 'cmake -DCMAKE_BUILD_TYPE=Debug'

cmake_minimum_required(VERSION 2.6)

project(bwtc)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS program_options)

if(PROFILER MATCHES 1)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -DPROFILER_ON")
else(PROFILER_ON MATCHES 1)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
endif()

if(ENTROPY_PROFILER MATCHES 1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENTROPY_PROFILER")
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions(-DNDEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -O3")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fno-inline")
endif()

set(EXECUTABLE_OUTPUT_PATH bin/)
set(OBJECT_FILE_PATH ${bwtc_SOURCE_DIR}/${EXECUTABLE_OUTPUT_PATH})

set(COMMON_SRC Coders.cpp BitCoders.cpp Utils.cpp Streams.cpp 
  WaveletCoders.cpp EntropyCoders.cpp HuffmanCoders.cpp PrecompressorBlock.cpp 
  BWTBlock.cpp)
add_library(common ${COMMON_SRC})

set(PREPROCESSORS_SRC_DIR preprocessors)
set(BWTRANSFORMS_SRC_DIR bwtransforms)
set(PROBMODELS_SRC_DIR probmodels)
set(TEST_SRC_DIR test)

add_subdirectory(${PREPROCESSORS_SRC_DIR})
add_subdirectory(${BWTRANSFORMS_SRC_DIR})
add_subdirectory(${PROBMODELS_SRC_DIR})

set(COMPRESSOR_SRC compress.cpp Compressor.cpp)
set(DECOMPRESSOR_SRC uncompress.cpp Decompressor.cpp)
set(PREPROCESSOR_SRC preprocess.cpp)
set(POSTPROCESSOR_SRC postprocess.cpp)

link_directories(${Boost_LIBRARY_DIR})
include_directories(${Boost_INCLUDE_DIR})

add_executable(compress ${COMPRESSOR_SRC})
target_link_libraries(compress common preprocessors probmodels bwtransforms 
  ${Boost_LIBRARIES})

add_executable(uncompress ${DECOMPRESSOR_SRC})
target_link_libraries(uncompress common preprocessors probmodels bwtransforms
  ${Boost_LIBRARIES})

add_executable(preprocess ${PREPROCESSOR_SRC})
target_link_libraries(preprocess common preprocessors 
  ${Boost_LIBRARIES})

add_executable(postprocess ${POSTPROCESSOR_SRC})
target_link_libraries(postprocess common preprocessors 
  ${Boost_LIBRARIES})

enable_testing()

set(CMAKE_BUILD_TYPE Debug)
add_subdirectory(${TEST_SRC_DIR})
