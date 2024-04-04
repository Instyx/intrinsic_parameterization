cmake_minimum_required(VERSION 3.16)
#set(CMAKE_CXX_STANDARD 17)
project(iparam)
set(CMAKE_PREFIX_PATH "/usr/inlcude/suitesparse" ${CMAKE_PREFIX_PATH})

# include(FetchContent)
# FetchContent_Declare(
#   tinyad
#   GIT_REPOSITORY https://github.com/patr-schm/tinyad.git
#   GIT_TAG 29417031c185b6dc27b6d4b684550d844459b735
# )
# #FetchContent_MakeAvailable(tinyad)
#
# FetchContent_GetProperties(tinyad)
# if(NOT tinyad_POPULATED)
#   # Fetch the content using previously declared details
#   FetchContent_Populate(tinyad)
#   message(STATUS "tinyad_SOURCE_DIR: ${tinyad_SOURCE_DIR}")
#   message(STATUS "tinyad_BINARY_DIR: ${tinyad_BINARY_DIR}")
#   add_subdirectory(${tinyad_SOURCE_DIR} ${tinyad_BINARY_DIR})
# endif()


list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_GLFW         "Build target igl::glfw"         ON)
option(LIBIGL_IMGUI        "Build target igl::imgui"        ON)
option(LIBIGL_OPENGL       "Build target igl::opengl"       ON)
#option(LIBIGL_PNG          "Build target igl::png"          ON)
#include(libigl)

add_subdirectory("libigl")
add_subdirectory("geometry-central")



# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
#igl_include(png)

# find_package (SuiteSparse NO_MODULE)
find_package(CHOLMOD REQUIRED)
# Add your project files
# file(GLOB SRC_FILES src/*.cpp)
file(GLOB IPARAM_FILES iparam/*.cpp)
file(GLOB IPARAM_HEADERS iparam/*.hpp)
add_library(iparam_core STATIC ${IPARAM_FILES} ${IPARAM_HEADERS})
add_library(predicates STATIC predicates.c)
target_include_directories(iparam_core PUBLIC iparam ${CHOLMOD_INCLUDES})
target_link_libraries(iparam_core PUBLIC igl::glfw igl::opengl igl::imgui geometry-central predicates SuiteSparse::CHOLMOD)

add_executable(cli src/cli.cpp ${IPARAM_HEADERS})
target_link_libraries(cli PUBLIC iparam_core)
target_compile_options(cli PRIVATE -O3)

add_executable(${PROJECT_NAME} src/main.cpp ${IPARAM_HEADERS})
target_link_libraries(${PROJECT_NAME} PUBLIC iparam_core)
target_compile_options(${PROJECT_NAME} PRIVATE -O3)
