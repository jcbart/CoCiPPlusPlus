cmake_minimum_required(VERSION 3.10)

project(CoCiP++ LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "--------------------------------------------")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "--------------------------------------------")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Submodules
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
add_subdirectory(submodules/yaml-cpp EXCLUDE_FROM_ALL)

# Target definition
file(GLOB SOURCES "src/*.cpp")
add_executable(CoCiP ${SOURCES})

target_include_directories(CoCiP
    PUBLIC include
)

target_link_libraries(CoCiP PRIVATE yaml-cpp)

# Compile features and options
target_compile_features(CoCiP PRIVATE cxx_std_11)
target_compile_options(CoCiP PRIVATE
    -Wall -Wextra
)
