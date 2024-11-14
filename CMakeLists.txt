cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(ALFI VERSION 0.0.0.0)

add_subdirectory(ALFI)

add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)