set(IPARAM_INCLUDE_DIR "../include")


# Define libraries
file(GLOB IPARAM_FILES *.cpp)
file(GLOB IPARAM_HEADERS ${IPARAM_INCLUDE_DIR}/*.hpp)

add_library(predicates STATIC predicates.c)

add_library(iparam_core STATIC ${IPARAM_FILES} ${IPARAM_HEADERS})
target_include_directories(iparam_core PUBLIC ${IPARAM_INCLUDE_DIR})
target_link_libraries(iparam_core PUBLIC igl::glfw igl::opengl igl::imgui geometry-central predicates)


