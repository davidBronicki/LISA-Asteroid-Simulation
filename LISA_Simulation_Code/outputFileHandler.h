#pragma once

#include <vector>
#include <string>

#include "vec3d.cuh"

namespace lisa{
	class outputFileHandler
	{
		std::vector<cudaUtil::vec3d> relativeAccelData;
		std::vector<std::vector<cudaUtil::vec3d>> absoluteAccelData;
		std::vector<std::vector<cudaUtil::vec3d>> unitVectorData;
	public:
		outputFileHandler();
		outputFileHandler(
			std::vector<cudaUtil::vec3d> inputRelData,
			std::vector<std::vector<cudaUtil::vec3d>> inputAbsData,
			std::vector<std::vector<cudaUtil::vec3d>> inputUnitVectors);
		void writeDataToFile(std::string fileDiscriminator,
			std::string headerString,
			std::vector<double> fileTimes,
			unsigned int seed);
	};
}
