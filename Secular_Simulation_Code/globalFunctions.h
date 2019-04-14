#pragma once

#include <vector>
#include <string>

#include "Orbit.h"

#include "defineOrbitalValues.h"

namespace lisa{
    float yearToEpoch(float inputYear);
    float epochToYear(float inputEpoch);
}

void writeCSVfile(std::vector<std::vector<double>> data, std::string location);

template<class T>
std::vector<T> concat(std::vector<T> a, std::vector<T> b){
    std::vector<T> output = std::vector<T>();
    for (int i = 0; i < a.size(); i++){
        output.push_back(a[i]);
    }
    for (int i = 0; i < b.size(); i++){
        output.push_back(b[i]);
    }
    return output;
}

inline double radians(double degrees)
{
    return degrees * PI / 180;
}

inline double degrees(double radians)
{
    return radians * 180 / PI;
}

#include "undefineOrbitalValues.h"
