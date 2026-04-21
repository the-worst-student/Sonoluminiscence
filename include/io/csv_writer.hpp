#pragma once

#include <string>
#include <vector>

<<<<<<< HEAD
#include "bubble/bubble_rhs.hpp"

class CsvWriter {
public:
    static void WriteBubbleSamples(const std::string& path, const std::vector<BubbleSample>& samples);
=======
#include <bubble/bubble_rhs.hpp>

class CsvWriter {
public:
    static void WriteBubbleSamples(
        const std::string& path,
        const std::vector<BubbleSample>& samples);
>>>>>>> 7de201c (Update project version)
};