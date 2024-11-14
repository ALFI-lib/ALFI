set(CMAKE_CXX_STANDARD 23)

find_package(GTest REQUIRED)

add_executable(test_dist dist/test_dist.cpp)

target_include_directories(test_dist PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_dist PRIVATE ALFI GTest::GTest GTest::Main)
target_compile_options(test_dist PRIVATE -Wall -Wextra -pedantic)

add_test(NAME test_dist COMMAND test_dist --gtest_color=yes)