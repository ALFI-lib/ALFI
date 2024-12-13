set(CMAKE_CXX_STANDARD 23)

find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(QCustomPlot REQUIRED)

add_executable(dist_QCustomPlot dist_QCustomPlot.cpp)
target_link_libraries(dist_QCustomPlot PRIVATE ALFI Qt5::Widgets QCustomPlot)
target_compile_options(dist_QCustomPlot PRIVATE -Wall -Wextra -pedantic -O3 -flto=auto -march=native)