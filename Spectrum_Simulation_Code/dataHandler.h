#pragma once

#include <vector>
#include <string>

#include "Orbit.cuh"
#include "LISA.h"
#include "specArgs.h"
#include "inputFileHandler.h"
#include "outputFileHandler.h"

namespace lisa{
	class dataHandler
	{
		specArgs argsHandler;
		inputFileHandler reader;
		outputFileHandler writer;
		std::vector<Orbit> simulationBodies;
		std::vector<double> generateFileTimes();
		std::string generateFileName();
		bool allIsWell;
	public:
		dataHandler();
		dataHandler(int argc, char** argv);
		std::vector<Orbit> generateAsteroids();
		LISA generateLISA();
		std::vector<float> generateSampleTimes();
		void handOffResults(
			std::vector<cudaUtil::vec3d> inputRelData,
			std::vector<std::vector<cudaUtil::vec3d>> inputAbsData,
			std::vector<std::vector<cudaUtil::vec3d>> inputUnitVectors,
			bool successfulSimulation);
		bool writeToFile();
		bool isGood();
		void killProgram();
	};
}
