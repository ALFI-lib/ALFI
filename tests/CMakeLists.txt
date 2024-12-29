set(CMAKE_CXX_STANDARD 23)

find_package(GTest REQUIRED)

macro(add_test_executable test_name test_source)
	add_executable(${test_name} ${test_source})
	target_include_directories(${test_name} PRIVATE ${GTEST_INCLUDE_DIRS})
	target_link_libraries(${test_name} PRIVATE ALFI GTest::GTest GTest::Main)
	target_compile_options(${test_name} PRIVATE -fsanitize=address,leak,undefined)
	target_link_options(${test_name} PRIVATE -fsanitize=address,leak,undefined)
	add_test(NAME ${test_name} COMMAND ${test_name} --gtest_color=yes)
endmacro()

add_test_executable(test_dist dist/test_dist.cpp)

add_test_executable(test_lagrange poly/test_lagrange.cpp)
add_test_executable(test_newton poly/test_newton.cpp)

add_test_executable(test_barycentric misc/test_barycentric.cpp)

add_test_executable(test_step spline/test_step.cpp)
add_test_executable(test_linear spline/test_linear.cpp)
add_test_executable(test_quadratic spline/test_quadratic.cpp)
add_test_executable(test_cubic spline/test_cubic.cpp)