cmake_minimum_required(VERSION 3.10)
project(solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# add_compile_options(-Wno-deprecated-declarations)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic -std=c++17")

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(Eigen3 REQUIRED NO_MODULE)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    $ENV{HOMEBREW_INCLUDE_PATH}
    $ENV{HOMEBREW_OPT_PATH}/eigen/include/eigen3
    $ENV{HOMEBREW_OPT_PATH}/boost/include
    $ENV{HOMEBREW_OPT_PATH}/ginac/include/ginac
    # $ENV{EXTERNAL_LIBS_PATH}/gnuplot/include # optional (currently not using)
)

message("current working directory ${CMAKE_CURRENT_SOURCE_DIR}")

link_directories(
    $ENV{HOMEBREW_OPT_PATH}/ginac/lib 
    $ENV{HOMEBREW_OPT_PATH}/boost/lib
)


file(GLOB_RECURSE SRC_FILES "${SRC_DIR}/*.cpp")

add_executable(main ${SRC_FILES} ${CMAKE_SOURCE_DIR}/main.cpp)
target_link_libraries(main PRIVATE ginac Boost::iostreams)

