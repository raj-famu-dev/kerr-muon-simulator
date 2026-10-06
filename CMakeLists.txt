cmake_minimum_required(VERSION 3.16)
project(Monograph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(monograph_core
        src/kerr.cpp
        src/rk4.cpp
        src/orbits.cpp
        src/csv_writer.cpp)
target_include_directories(monograph_core PUBLIC include)

add_executable(monograph_isco   src/main_isco.cpp)
add_executable(monograph_plunge src/main_plunge.cpp)
target_link_libraries(monograph_isco   PRIVATE monograph_core)
target_link_libraries(monograph_plunge PRIVATE monograph_core)

enable_testing()
foreach(t test_kerr test_rk4 test_isco test_constraint)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE monograph_core)
    add_test(NAME ${t} COMMAND ${t})
endforeach()