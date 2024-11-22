cmake_minimum_required(VERSION 3.14)

# Project name and version
project(MDExample VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED)

# Add the executable
add_executable(md src/main.cpp)

target_link_libraries(md dl)
