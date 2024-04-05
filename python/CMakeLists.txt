file(GLOB PYTHON_FILES python/*.cpp)

nanobind_add_module(pyiparam ${PYTHON_FILES} module.cpp)
target_link_libraries(pyiparam PUBLIC iparam_core)
