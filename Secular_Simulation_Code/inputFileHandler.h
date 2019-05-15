#pragma once

#include <vector>
#include <string>

#include "Orbit.h"

namespace lisa{
	class inputFileHandler{
		int currentRandIndex;
		unsigned int seed;
		std::vector<float> orbitalData;
		std::vector<std::string> orbitNames;
		std::vector<float> normalRands;
		std::vector<float> getOrbitData(int index);
		float getRandValue(float mean, float standardDeviation);
		float findMass(float absMag, float type, float albedo);
	public:
		inputFileHandler();
		inputFileHandler(int numberOfRandomNumbers);
		Orbit generateOrbit(int index, float time0);
		std::vector<Orbit> generateOrbits(std::vector<int> indices, float time0);
		std::string getName(int index);
		unsigned int getSeed();
	};
}