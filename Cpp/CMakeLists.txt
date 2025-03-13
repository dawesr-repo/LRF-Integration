cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_COMPILER "g++") 
set(CMAKE_MAKE_PROGRAM "make")

project(lrf_integration LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)

add_executable(lrf_integration
    main.cpp
    evaluate_lrf.cpp
    interactions.cpp
    calculate_tensor.cpp
    potential_energy_surface.cpp
    potential_energy_surface.h    
)
 