#pragma once

#include <vector>
#include <string>

#include "argumentsHandler.h"

namespace lisa{
	class specArgs
	{
		arguments::argumentsHandler handler;
		// int randCount, range0, range1, sampleTimeCount,
		// 	dataLength;
		// float dt, t0, t1;
	public:
		specArgs();
		specArgs(int argc, char** argv);
		std::vector<int> getSelection();
		std::vector<float> generateSampleTimes();
		float getInitialEpoch();
		std::string getFileName(unsigned int seed);
		int getRandCount();
		std::string generateHeaderString();
		std::string getComment();
		bool isGood();
	};
}
