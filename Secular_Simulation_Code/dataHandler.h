#pragma once

#include <vector>
#include <string>

#include "Orbit.h"
#include "LISA.h"
#include "specArgs.h"
#include "inputFileHandler.h"

namespace lisa{
	class dataHandler
	{
		specArgs argsHandler;
		inputFileHandler reader;
		std::vector<Orbit> simulationBodies;
		bool allIsWell;
	public:
		dataHandler();
		dataHandler(int argc, char** argv);
		std::vector<Orbit> generateAsteroids();
		LISA generateLISA();
		bool isGood();
	};
}
