#pragma once

#include <chrono>

#define ADD_TIME_BLOCK(times, ...) do {                                           \
    auto __start = std::chrono::high_resolution_clock::now();                    \
    __VA_ARGS__                                                                  \
    auto __end = std::chrono::high_resolution_clock::now();                      \
	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
	times += __elapsed.count();    \
} while(0);

#define ADD_TO_VECTOR_TIME_BLOCK(times, ...) do {                                           \
    auto __start = std::chrono::high_resolution_clock::now();                    \
    __VA_ARGS__                                                                  \
    auto __end = std::chrono::high_resolution_clock::now();                      \
	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
	times.push_back(__elapsed.count());    \
} while(0);
// #define ADD_TIME_BLOCK(times, ...) do {                                           \
//     auto __start = std::chrono::high_resolution_clock::now();                    \
//     __VA_ARGS__                                                                  \
//     auto __end = std::chrono::high_resolution_clock::now();                      \
// 	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
// 	times.add(__elapsed.count());    \
// } while(0);

#define TIME_BLOCK(name, ...) do {                                           \
    auto __start = std::chrono::high_resolution_clock::now();                    \
    __VA_ARGS__                                                                  \
    auto __end = std::chrono::high_resolution_clock::now();                      \
	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
	std::cout << name << " took: " << __elapsed.count() << "ms" << std::endl;    \
} while(0);
