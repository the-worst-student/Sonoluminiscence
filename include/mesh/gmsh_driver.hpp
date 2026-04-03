#pragma once

#include <string>

#include "core/config.hpp"

class GmshDriver {
public:
    explicit GmshDriver(const ProjectConfig& config);

    void BuildAxisymmetricMesh(const std::string& output_mesh_path) const;

private:
    ProjectConfig config_;
};