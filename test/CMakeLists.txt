cmake_minimum_required(VERSION 2.6)

project(bwtc_tests)

set(EXECUTABLE_OUTPUT_PATH bin/)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS system filesystem)

link_directories(${Boost_LIBRARY_DIR} ${OBJECT_FILE_PATH})
include_directories(${Boost_INCLUDE_DIR})

set(BOOST_UNIT_TESTS UtilsTest WaveletTest PairReplacerTest GrammarTest)

foreach(program ${BOOST_UNIT_TESTS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} common boost_unit_test_framework preprocessors)
  add_test(${program} ${EXECUTABLE_OUTPUT_PATH}/${program})
  set_tests_properties(${program} PROPERTIES FAIL_REGULAR_EXPRESSION 
    "[.\n]*failure")
endforeach(program)

add_executable(CompressorAndDecompressorTest CompressorAndDecompressorTest.cpp
  ../Decompressor.cpp ../Compressor.cpp)
target_link_libraries(CompressorAndDecompressorTest 
  common boost_unit_test_framework preprocessors bwtransforms probmodels)
add_test(CompressorAndDecompressorTest 
  ${EXECUTABLE_OUTPUT_PATH}/CompressorAndDecompressorTest)
set_tests_properties(CompressorAndDecompressorTest 
  PROPERTIES FAIL_REGULAR_EXPRESSION "[.\n]*failure")

add_executable(SaisTest SaisTest.cpp)
target_link_libraries(SaisTest common)
add_test(SaisTest ${EXECUTABLE_OUTPUT_PATH}/SaisTest)
set_tests_properties(SaisTest PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(DivsufsortTest DivsufsortTest.cpp)
target_link_libraries(DivsufsortTest bwtransforms)

add_executable(RawStreamTest RawStreamTest.cpp)
target_link_libraries(RawStreamTest common ${Boost_LIBRARIES})
add_test(RawStreamTest ${EXECUTABLE_OUTPUT_PATH}/RawStreamTest "test_file")
set_tests_properties(RawStreamTest PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(PrecompressorTest PrecompressorTest.cpp)
target_link_libraries(PrecompressorTest common preprocessors)
add_test(PrecompressorTest ${EXECUTABLE_OUTPUT_PATH}/PrecompressorTest "test_file")
set_tests_properties(PrecompressorTest PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(InverseBwtTest InverseBwtTest.cpp)
target_link_libraries(InverseBwtTest common bwtransforms)

add_executable(InverseBwtOnFileTest InverseBwtOnFileTest.cpp)
target_link_libraries(InverseBwtOnFileTest common bwtransforms)

add_executable(LFpowersTest LFpowersTest.cpp)
target_link_libraries(LFpowersTest common bwtransforms)
