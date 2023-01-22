cmake_minimum_required(VERSION 3.7)
project(algorithm)

set(CMAKE_CXX_STANDARD 23)

add_executable(algorithm main.cpp main.h Node.cpp Node.h LabPixel.cpp LabPixel.h CLIInputParser.h)

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
target_link_libraries( algorithm ${OpenCV_LIBS} )

find_package(SDL2 REQUIRED COMPONENTS SDL2)
target_link_libraries(algorithm SDL2::SDL2)


