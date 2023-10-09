cmake_minimum_required(VERSION 3.16)
project(example)




list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_GLFW         "Build target igl::glfw"         ON)
option(LIBIGL_IMGUI        "Build target igl::imgui"        ON)
option(LIBIGL_OPENGL       "Build target igl::opengl"       ON)
option(LIBIGL_PNG          "Build target igl::png"          ON)
include(libigl)

add_subdirectory("geometry-central")




# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
igl_include(png)


# Add your project files
file(GLOB SRC_FILES src/*.cpp)
file(GLOB IPARAM_FILES iparam/*.cpp)
file(GLOB IPARAM_HEADERS iparam/*.hpp)
add_library(iparam_core STATIC ${IPARAM_FILES} ${IPARAM_HEADERS})
target_include_directories(iparam_core PUBLIC iparam)
target_link_libraries(iparam_core PUBLIC igl::glfw igl::opengl igl::imgui igl::png geometry-central)

add_executable(${PROJECT_NAME} ${SRC_FILES} ${IPARAM_HEADERS})
target_link_libraries(${PROJECT_NAME} PUBLIC iparam_core) 
target_compile_options(${PROJECT_NAME} PRIVATE -O3 ) 
