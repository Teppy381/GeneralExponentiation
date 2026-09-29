find_package(doctest CONFIG QUIET)
if(NOT doctest_FOUND)
    message(STATUS "doctest not found: tests are skipped")
    return()
endif()

# static_checks.cpp holds only static_asserts: building gp_tests is half of the test.
add_executable(gp_tests static_checks.cpp runtime_tests.cpp)
target_link_libraries(gp_tests PRIVATE gp gp_options doctest::doctest)
# doctest 2.4 vs clang 22 / libstdc++ 16: __COUNTER__ in its macros and #include <ciso646>.
# Source-level options come after -Wpedantic from gp_options, so they stick.
set_source_files_properties(runtime_tests.cpp PROPERTIES COMPILE_OPTIONS "-Wno-c2y-extensions;-Wno-#warnings")
add_test(NAME gp_tests COMMAND gp_tests)
