cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 20)

project(phisics
        VERSION 0.1
        LANGUAGES CXX
        )

include(cmake/conan.cmake)
run_conan()

add_executable(main main2.cpp include/lab4.hpp include/common.hpp include/lab2.hpp include/lab3.hpp)

target_link_libraries(main PRIVATE ${CONAN_LIBS})

