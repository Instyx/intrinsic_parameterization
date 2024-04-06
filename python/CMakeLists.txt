file(GLOB PYTHON_FILES *.cpp)

nanobind_add_module(pyiparam ${PYTHON_FILES} module.cpp)
target_link_libraries(pyiparam PUBLIC iparam_core)
