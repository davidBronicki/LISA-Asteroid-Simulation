#pragma once

#include <vector>
#include <string>

#include "Orbit.cuh"
#include "vec3d.cuh"

namespace lisa{
    float yearToEpoch(float inputYear);
    float epochToYear(float inputEpoch);
    std::vector<double> simTimeToFileTime(std::vector<float> simTime, float initialEpoch);
    std::string generateHeaderString(
        unsigned int seed,
        std::vector<int> asteroidIndices,
        std::vector<lisa::Orbit> bodiesSimulated,
        std::vector<float> timesSampled,
        float initialYear,
        std::string comment);
    std::vector<double> floatVectorToDoubleVector(std::vector<float> input);
    std::vector<double> vecVectorToDoubleVector(std::vector<cudaUtil::vec3d> input);
}

void writeCSVfile(std::vector<std::vector<double>> data, std::string location);
std::string str(double input);

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