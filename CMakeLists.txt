cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 20)

project(phisics
        VERSION 0.1
        LANGUAGES CXX
        )

include(cmake/conan.cmake)
run_conan()

add_executable(main main.cpp)

target_link_libraries(main PRIVATE ${CONAN_LIBS})

