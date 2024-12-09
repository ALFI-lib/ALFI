set(CMAKE_CXX_STANDARD 23)

find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(QCustomPlot REQUIRED)

add_executable(interpolation interpolation.cpp)
target_link_libraries(interpolation PRIVATE ALFI Qt5::Widgets QCustomPlot)
target_compile_options(interpolation PRIVATE -Wall -Wextra -pedantic -O3 -flto -march=native)