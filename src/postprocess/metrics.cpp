#pragma once

#include <string>
#include <vector>

#include "bubble/bubble_rhs.hpp"

class CsvWriter {
 public:
    static void WriteBubbleSamples(const std::string& path, const std::vector<BubbleSample>& samples);
};
