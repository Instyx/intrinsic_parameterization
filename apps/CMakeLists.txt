add_executable(cli cli.cpp)
target_link_libraries(cli PUBLIC iparam_core)
target_compile_options(cli PRIVATE -O3)

add_executable(main main.cpp)
target_link_libraries(main PUBLIC iparam_core)
target_compile_options(main PRIVATE -O3)

add_executable(filter filter.cpp)
target_link_libraries(filter PUBLIC iparam_core)
target_compile_options(filter PRIVATE -O3)
