cmake_minimum_required(VERSION 2.6)

project(probmodels)

FILE(GLOB cppSourceFiles *.cpp *.cc)
FILE(GLOB hppHeaders *.hpp *.h)

set(PROBMODEL_SOURCES ${cppSourceFiles} ${hppHeaders})

add_library(probmodels ${cppSourceFiles})
