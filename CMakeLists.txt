cmake_minimum_required(VERSION 3.16)

project(LRF_Integration LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_executable(LRF_Integration
                    Cpp/main.cpp
                    Cpp/PotentialEnergySurface.h
                    Cpp/PotentialEnergySurface.cpp
                    Cpp/interactions.cpp
        Cpp/tensor.cpp
        Cpp/tensor.h
)
target_link_libraries(LRF_Integration)


#How to run it
#cmake -S . -B build
#cmake --build build -j 12
# for testing only: cd build && ctest --output-on-failure -j12