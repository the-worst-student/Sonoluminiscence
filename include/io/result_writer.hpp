#pragma once

#include <string>

class ResultWriter {
public:
    static void WriteTextFile(
        const std::string& output_path,
        const std::string& content);
};
