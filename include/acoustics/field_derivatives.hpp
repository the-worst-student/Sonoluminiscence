#pragma once

#include <acoustics/field_sampler.hpp>

class FieldDerivatives {
public:
    static void AddPressureGradient(AcousticFieldData* field_data);
};
