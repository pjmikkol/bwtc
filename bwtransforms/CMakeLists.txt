cmake_minimum_required(VERSION 2.6)

project(bwtranforms)

FILE(GLOB cppSourceFiles *.cpp *.cc)
FILE(GLOB hppHeaders *.hpp *.h)

set(BWT_SOURCES ${cppSourceFiles} ${hppHeaders})

add_library(bwtransforms ${cppSourceFiles})
