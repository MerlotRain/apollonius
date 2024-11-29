cmake_minimum_required(VERSION 3.10)
project(apollonius LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 99)

option(BUILD_TESTING "build apollonius solution tests" OFF)

add_library(apollonius STATIC apollonius.cpp)

if(BUILD_TESTING)

file(GLOB test_files ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
foreach(file ${test_files})
    
endforeach()
endif()
