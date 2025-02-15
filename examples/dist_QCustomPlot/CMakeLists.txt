set(CMAKE_CXX_STANDARD 17)

find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(QCustomPlot REQUIRED)

add_executable(dist_QCustomPlot dist_QCustomPlot.cpp)
target_link_libraries(dist_QCustomPlot PRIVATE ALFI Qt5::Widgets QCustomPlot)