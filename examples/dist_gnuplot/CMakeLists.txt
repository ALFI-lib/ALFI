set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(Gnuplot REQUIRED)
find_path(GNUPLOT_IOSTREAM_INCLUDE gnuplot-iostream.h REQUIRED)

add_executable(dist_gnuplot dist_gnuplot.cpp)
target_include_directories(dist_gnuplot PRIVATE ${Boost_INCLUDE_DIRS} ${GNUPLOT_IOSTREAM_INCLUDE})
target_link_libraries(dist_gnuplot PRIVATE ALFI ${Boost_LIBRARIES})