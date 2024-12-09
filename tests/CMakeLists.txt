set(CMAKE_CXX_STANDARD 23)

find_package(GTest REQUIRED)

add_executable(test_dist dist/test_dist.cpp)
target_include_directories(test_dist PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_dist PRIVATE ALFI GTest::GTest GTest::Main)
target_compile_options(test_dist PRIVATE -Wall -Wextra -pedantic)
add_test(NAME test_dist COMMAND test_dist --gtest_color=yes)

add_executable(test_lagrange poly/test_lagrange.cpp)
target_include_directories(test_lagrange PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_lagrange PRIVATE ALFI GTest::GTest GTest::Main)
target_compile_options(test_lagrange PRIVATE -Wall -Wextra -pedantic)
add_test(NAME test_lagrange COMMAND test_lagrange --gtest_color=yes)

add_executable(test_newton poly/test_newton.cpp)
target_include_directories(test_newton PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_newton PRIVATE ALFI GTest::GTest GTest::Main)
target_compile_options(test_newton PRIVATE -Wall -Wextra -pedantic)
add_test(NAME test_newton COMMAND test_newton --gtest_color=yes)

add_executable(test_barycentric bary/test_barycentric.cpp)
target_include_directories(test_barycentric PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(test_barycentric PRIVATE ALFI GTest::GTest GTest::Main)
target_compile_options(test_barycentric PRIVATE -Wall -Wextra -pedantic)
add_test(NAME test_barycentric COMMAND test_barycentric --gtest_color=yes)