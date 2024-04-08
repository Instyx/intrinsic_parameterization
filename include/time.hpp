#pragma once

#include <array>
#include <cstddef>

/**
 * Ring buffer for recording timings
 */
class Times
{
public:
	Times();

	void add(float time);

	float const* data() const;

	size_t size() const;

private:
	std::array<float, 32u> m_times;
};

struct Timings
{
	Times precompute;
	Times subdivide;
};

#define ADD_TIME_BLOCK(times, ...) do {                                           \
    auto __start = std::chrono::high_resolution_clock::now();                    \
    __VA_ARGS__                                                                  \
    auto __end = std::chrono::high_resolution_clock::now();                      \
	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
	times.add(__elapsed.count());    \
} while(0);

#define TIME_BLOCK(name, ...) do {                                           \
    auto __start = std::chrono::high_resolution_clock::now();                    \
    __VA_ARGS__                                                                  \
    auto __end = std::chrono::high_resolution_clock::now();                      \
	std::chrono::duration<float, std::milli> __elapsed = __end - __start; \
	std::cout << name << " took: " << __elapsed.count() << "ms" << std::endl;    \
} while(0);
