cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 20)

project(phisics)

add_executable(main main2.cpp)

target_include_directories(main PRIVATE 
	boost/libs/math/include
	boost/libs/multiprecision/include
	)

