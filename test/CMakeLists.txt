cmake_minimum_required(VERSION 2.6)

project(bwtc_tests)

set(EXECUTABLE_OUTPUT_PATH bin/)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS system)

link_directories(${Boost_LIBRARY_DIR} ${OBJECT_FILE_PATH})
include_directories(${Boost_INCLUDE_DIR})

set(BOOST_UNIT_TESTS utils_test expander_test prime_sequence_detector_test
  mask_sequence_detector_test)

foreach(program ${BOOST_UNIT_TESTS})
  add_executable(${program} ${program}.cc)
  target_link_libraries(${program} common boost_unit_test_framework preprocessors)
  add_test(${program} ${CMAKE_CURRENT_BINARY_DIR}/${program})
  set_tests_properties(${program} PROPERTIES FAIL_REGULAR_EXPRESSION 
    "[.\n]*failure")
endforeach(program)

add_executable(coders_test coders_test.cc)
target_link_libraries(coders_test common probmodels)
add_test(coders_test ${EXECUTABLE_OUTPUT_PATH}/coders_test "test_file")
set_tests_properties(coders_test PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

# TODO: Poista huono toteutus
add_executable(sa-is_test sa-is_test.cc)
target_link_libraries(sa-is_test common)
add_test(sa-is_test ${EXECUTABLE_OUTPUT_PATH}/sa-is_test)
set_tests_properties(sa-is_test PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(stream_test stream_test.cc)
target_link_libraries(stream_test common ${Boost_LIBRARIES})
add_test(stream_test ${EXECUTABLE_OUTPUT_PATH}/stream_test "test_file")
set_tests_properties(stream_test PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(preproc_test preproc_test.cc)
target_link_libraries(preproc_test common preprocessors)
add_test(preproc_test ${EXECUTABLE_OUTPUT_PATH}/preproc_test "test_file")
set_tests_properties(preproc_test PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")

add_executable(preproc_algo_test preproc_algo_test.cc)
target_link_libraries(preproc_algo_test common preprocessors)

add_executable(bwt_and_preproctest bwt_and_preproctest.cc)
target_link_libraries(bwt_and_preproctest common preprocessors bwtransforms)

add_executable(long_sequences_test long_sequences_test.cc)
target_link_libraries(long_sequences_test common preprocessors)



#add_test(long_sequences_test ${EXECUTABLE_OUTPUT_PATH}/long_sequences_test)
#set_tests_properties(long_sequences_test PROPERTIES PASS_REGULAR_EXPRESSION ".*pass")