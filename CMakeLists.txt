cmake_minimum_required(VERSION 3.10)

project(CoCiP++ LANGUAGES CXX)

# Build type

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

option(COCIP_BUILD_EXECUTABLE "Build a standalone executable" ON)

message(STATUS "--------------------------------------------")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Build a standalone executable: ${COCIP_BUILD_EXECUTABLE}")
message(STATUS "--------------------------------------------")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Submodules

if(NOT TARGET yaml-cpp)
    set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
    add_subdirectory(submodules/yaml-cpp EXCLUDE_FROM_ALL)
endif()

# Compilation config interface
# (specifies settings for both library and executable)

add_library(build_config INTERFACE)
target_compile_features(build_config INTERFACE cxx_std_11)
target_compile_options(build_config INTERFACE -Wall -Wextra)
target_link_libraries(build_config INTERFACE
    yaml-cpp
)

# Library sources

set(LIB_SOURCES
    "src/CoCiP.cpp"
    "src/CoCiPLog.cpp"
    "src/contrail_properties.cpp"
    "src/met.cpp"
    "src/params.cpp"
    "src/radiative_forcing.cpp"
    "src/radiative_heating.cpp"
    "src/sac.cpp"
    "src/wake_vortex.cpp"
)

# Library definition

add_library(cocip STATIC ${LIB_SOURCES})
target_include_directories(cocip PUBLIC include)
target_link_libraries(cocip PRIVATE build_config)

# Executable definition

if(COCIP_BUILD_EXECUTABLE)
    add_executable(CoCiP "src/main.cpp")
    target_link_libraries(CoCiP PRIVATE build_config cocip)
endif()

# Copy default input files if they do not already exist in build directory
set(INPUT_FILES "CoCiP-config.yaml" "CoCiP-params.yaml")

foreach(FILE_NAME ${INPUT_FILES})
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}")
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/inputs/${FILE_NAME}.default"
            "${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}"
            COPYONLY
        )
    endif()
endforeach()