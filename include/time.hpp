#pragma once

#include <chrono>
#include <utility>

namespace chrono = std::chrono;

auto funcTime =
    [](auto&& func, auto&&... params) {
        // get time before function invocation
        const auto& start = std::chrono::steady_clock::now();
        // function invocation using perfect forwarding
        std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
        // get time after function invocation
        const auto& stop = std::chrono::steady_clock::now();
        return chrono::duration_cast<chrono::microseconds>(stop - start).count();
     };
