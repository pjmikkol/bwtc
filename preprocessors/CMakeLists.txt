cmake_minimum_required(VERSION 2.6)

project(preprocessors)

FILE(GLOB cppSourceFiles *.cpp *.cc)
FILE(GLOB hppHeaders *.hpp *.h)

set(PREPROCESSOR_SOURCES ${cppSourceFiles} ${hppHeaders})

add_library(preprocessors ${cppSourceFiles})

target_link_libraries(preprocessors common)