#pragma once

#include <string>

#include <core/config.hpp>

class GeometryBuilder {
    public:
     explicit GeometryBuilder(const ProjectConfig& config);

     void BuildAxisymmetricGeometry() const;

    private:
     ProjectConfig config_;
};