set(CMAKE_CXX_STANDARD 23)

find_package(GTest REQUIRED)

set(DEPS_DIR "${CMAKE_CURRENT_LIST_DIR}/deps")
set(TEST_DATA_DIR "${CMAKE_CURRENT_LIST_DIR}/test_data")

macro(add_test_executable test_name)
	add_executable(${test_name} ${ARGN} test_utils.h)
	target_include_directories(${test_name} PRIVATE ${GTEST_INCLUDE_DIRS} ${DEPS_DIR})
	target_link_libraries(${test_name} PRIVATE ALFI GTest::GTest GTest::Main)
	target_compile_options(${test_name} PRIVATE -fsanitize=address,leak,undefined)
	target_link_options(${test_name} PRIVATE -fsanitize=address,leak,undefined)
	target_compile_definitions(${test_name} PRIVATE TEST_DATA_DIR="${TEST_DATA_DIR}")
	add_test(NAME ${test_name} COMMAND ${test_name} --gtest_color=yes)
endmacro()

add_test_executable(test_dist dist/test_dist.cpp)

add_test_executable(test_poly poly/test_poly.cpp)

add_test_executable(test_barycentric misc/test_barycentric.cpp)

add_test_executable(test_polyeqv spline/test_polyeqv.cpp)

add_test_executable(test_step spline/test_step.cpp)
add_test_executable(test_linear spline/test_linear.cpp)
add_test_executable(test_quadratic spline/test_quadratic.cpp)
add_test_executable(test_cubic spline/test_cubic.cpp)