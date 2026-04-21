#pragma once

#include <string>
#include <vector>

#include <bubble/bubble_rhs.hpp>

class BubbleVtkWriter {
public:
    struct Options {
        std::string OutputDirectory = "bubble_frames";
        int AzimuthalSegments = 48;
        int PolarSegments = 24;
    };

    static void WriteAnimationFrames(
        const std::vector<BubbleSample>& samples,
        double bubble_r_m,
        double bubble_z_m,
        const Options& options);
};
