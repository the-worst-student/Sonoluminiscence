#pragma once

enum class BoundaryTag : int {
    cWall = 1, // without roof, because source is on roof
    cReflector = 2,
    cSource = 3,
    cAxis = 4
};

enum class SurfaceTag : int {
    cFluid = 1,
    cReflectorBody = 2
};
