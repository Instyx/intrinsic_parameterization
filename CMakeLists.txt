cmake_minimum_required(VERSION 3.16..3.21)

project(iparam)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# Configure general CMake variables
set(CMAKE_DEBUG_POSTFIX d)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(IPARAM_EXTERNAL_DIR "external")


if (WIN32)
	# Enable parallel builds for Visual Studio and C++
	add_compile_options("$<$<AND:$<CXX_COMPILER_ID:MSVC>,$<COMPILE_LANGUAGE:CXX>>:/MP>")

	# Add central defines for windows
	add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES)
endif()

# Add Geometry Central
add_subdirectory(${IPARAM_EXTERNAL_DIR}/geometry-central)

# Add libigl
option(LIBIGL_GLFW         "Build target igl::glfw"         ON)
option(LIBIGL_IMGUI        "Build target igl::imgui"        ON)
option(LIBIGL_OPENGL       "Build target igl::opengl"       ON)
#option(LIBIGL_PNG          "Build target igl::png"          ON)

add_subdirectory(${IPARAM_EXTERNAL_DIR}/libigl)

igl_include(glfw)
igl_include(imgui)
igl_include(copyleft cgal)

# Add Cholmod
find_package(SuiteSparse REQUIRED)
include_directories(${CHOLMOD_INCLUDE_DIR})


# Add Nanobind
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind)

FetchContent_GetProperties(nanobind)
if(NOT nanobind)
  FetchContent_Populate(nanobind)
endif()

add_subdirectory(${nanobind_SOURCE_DIR})



#######################
# Our libraries and applications
#######################

# The compiled library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(apps)

# The executable python binding code is here
add_subdirectory(python)
