#pragma once

#include <string>

#include "core/config.hpp"

class YamlReader {
 public:
  static ProjectConfig ReadProjectConfig(const std::string& file_path);
};
