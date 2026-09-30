find_package(doctest CONFIG QUIET)
if(NOT doctest_FOUND)
    message(STATUS "doctest not found: tests are skipped")
    return()
endif()

# static_checks.cpp holds only static_asserts: building gp_tests is half of the test.
add_executable(gp_tests static_checks.cpp runtime_tests.cpp)
target_link_libraries(gp_tests PRIVATE gp doctest::doctest)
add_test(NAME gp_tests COMMAND gp_tests)
