cmake_minimum_required(VERSION 2.6)

project(bwtranforms)

FILE(GLOB cppSourceFiles *.cpp *.cc)
FILE(GLOB hppHeaders *.hpp *.h)

set(sources ${cppSourceFiles} ${hppHeaders})
set (EXECUTABLE_OUTPUT_PATH bin/)

add_library(bwtransforms ${cppSourceFiles})
