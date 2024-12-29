cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g3")
set(CMAKE_CXX_FLAGS_RELEASE "-g0 -s -O3 -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -O2 -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-g0 -s -Oz -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_FAST "-g0 -s -Ofast -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_SANITIZE "-g3 -fsanitize=address,leak,undefined")

project(ALFI VERSION 0.0.0.0)

add_subdirectory(ALFI)

add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)