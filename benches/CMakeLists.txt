set(CMAKE_CXX_STANDARD 23)

include(ExternalProject)

set(BENCHMARK_BUILD_DIR "${CMAKE_BINARY_DIR}/benches/deps/benchmark")

ExternalProject_Add(benchmark
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/benchmark"
	BINARY_DIR "${BENCHMARK_BUILD_DIR}"
	CMAKE_ARGS
		-DCMAKE_BUILD_TYPE=Release
		"-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
		-DBENCHMARK_ENABLE_INSTALL=OFF
		-DBENCHMARK_ENABLE_TESTING=OFF
		-DBENCHMARK_ENABLE_LTO=ON
	INSTALL_COMMAND ""
)

add_library(benchmark::benchmark STATIC IMPORTED)

set_target_properties(benchmark::benchmark PROPERTIES
	IMPORTED_LOCATION
		"${BENCHMARK_BUILD_DIR}/src/libbenchmark.a"
	INTERFACE_INCLUDE_DIRECTORIES
		"${CMAKE_CURRENT_SOURCE_DIR}/deps/benchmark/include"
)

add_dependencies(benchmark::benchmark benchmark)

macro(add_bench_executable bench_name)
	add_executable(${bench_name} ${ARGN} bench_utils.h bench_data.h)
	target_link_libraries(${bench_name} PRIVATE ALFI benchmark::benchmark)
endmacro()

add_bench_executable(bench_poly poly/bench_poly.cpp)
add_bench_executable(bench_ratf ratf/bench_ratf.cpp)
add_bench_executable(bench_barycentric misc/bench_barycentric.cpp)
add_bench_executable(bench_polyeqv spline/bench_polyeqv.cpp)
add_bench_executable(bench_step spline/bench_step.cpp)
add_bench_executable(bench_linear spline/bench_linear.cpp)
add_bench_executable(bench_quadratic spline/bench_quadratic.cpp)
add_bench_executable(bench_cubic spline/bench_cubic.cpp)