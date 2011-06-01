# To compile release version run  'cmake -DCMAKE_BUILD_TYPE=Release'
# For debug version run 'cmake -DCMAKE_BUILD_TYPE=Debug'

cmake_minimum_required(VERSION 2.6)

project(bwtc)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS program_options)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions(-DNDEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -O3")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()


set(EXECUTABLE_OUTPUT_PATH bin/)

set(PREPROCESSORS_SRC_DIR preprocessors)
set(BWTRANSFORMS_SRC_DIR bwtransforms)
set(PROBMODELS_SRC_DIR probmodels)

add_subdirectory(${PREPROCESSORS_SRC_DIR})
add_subdirectory(${BWTRANSFORMS_SRC_DIR})
add_subdirectory(${PROBMODELS_SRC_DIR})

# Fix coders.cpp
set(COMMON_SRC MainBlock.cpp BlockManager.cpp coders.cc rl_compress.cc utils.cc stream.cc)
set(COMPRESSOR_SRC compr.cc)
set(DECOMPRESSOR_SRC uncompr.cc)

link_directories(${Boost_LIBRARY_DIR})
include_directories(${Boost_INCLUDE_DIR})

add_executable(compr ${COMMON_SRC} ${COMPRESSOR_SRC})
target_link_libraries(compr preprocessors)
target_link_libraries(compr probmodels)
target_link_libraries(compr bwtransforms)
target_link_libraries(compr ${Boost_LIBRARIES})

add_executable(uncompr ${COMMON_SRC} ${DECOMPRESSOR_SRC})
target_link_libraries(uncompr preprocessors)
target_link_libraries(uncompr probmodels)
target_link_libraries(uncompr bwtransforms)
target_link_libraries(uncompr ${Boost_LIBRARIES})
