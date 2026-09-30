cmake_minimum_required(VERSION 3.28)
project(GeneralExponentiation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GP_BUILD_TESTS "Build the tests" ON)

# Eigen: the system package if there is one, otherwise the headers from GitLab.
find_package(Eigen3 CONFIG QUIET)
if(NOT TARGET Eigen3::Eigen)
    include(FetchContent)
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.1
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR no-cmake)  # headers only: skip Eigen's own build
    FetchContent_MakeAvailable(eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
endif()

# The header-only library
add_library(gp INTERFACE)
target_include_directories(gp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(gp INTERFACE cxx_std_23)
target_link_libraries(gp INTERFACE Eigen3::Eigen)

add_executable(showcase main.cpp)
target_link_libraries(showcase PRIVATE gp)

if(GP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
