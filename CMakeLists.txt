cmake_minimum_required(VERSION 3.10)
project(apollonius LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 99)

option(BUILD_TESTING "build apollonius solution tests" OFF)

add_library(apollonius STATIC apollonius.cpp)

if(BUILD_TESTING)
    add_executable(apollonius_test apollonius_test.cpp)
    target_link_libraries(apollonius_test apollonius)
endif()
