cmake_minimum_required(VERSION 3.16)
project(raytracer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -pthread -flto")
add_executable(raytracer raytracer.cpp raytracer.h parser.h parser.cpp tinyxml2.h tinysml2.cpp util.h util.cpp ppm.h ppm.cpp definitions.h)
