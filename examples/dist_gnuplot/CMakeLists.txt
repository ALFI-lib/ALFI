set(CMAKE_CXX_STANDARD 23)

find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(Gnuplot REQUIRED)

add_executable(dist_gnuplot dist_gnuplot.cpp)
target_include_directories(dist_gnuplot PRIVATE ${Boost_INCLUDE_DIRS} ${Gnuplot_INCLUDE_DIRS})
target_link_libraries(dist_gnuplot PRIVATE ALFI ${Boost_LIBRARIES} ${Gnuplot_LIBRARIES})
target_compile_options(dist_gnuplot PRIVATE -Wall -Wextra -pedantic -O3 -flto=auto -march=native)