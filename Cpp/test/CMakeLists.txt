cmake_minimum_required(VERSION 3.16)
project(TestUnit)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable testing for this project
include(CTest)

enable_testing()

add_executable(
        TestUnit
        testRunner.cpp
        ../sample2.cpp
        ../sample2.h
)
target_link_libraries(
        TestUnit
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(TestUnit)

#How to run it
# go to Cpp\test folder and run:
# cmake -S . -B build
# cmake --build build -j 12
# cd build && ctest --output-on-failure -j12