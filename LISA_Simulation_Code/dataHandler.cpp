#include "globalSwitchStatements.h"
#include "dataHandler.h"

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>

#include "Orbit.cuh"
#include "LISA.h"
#include "specArgs.h"
#include "inputFileHandler.h"
#include "outputFileHandler.h"
#include "globalFunctions.h"

#include "defineOrbitalValues.h"

using namespace std;
using namespace lisa;
using namespace cudaUtil;

#ifdef DEBUG_MODE
#define print(input) cout << "data: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

//The general data handler for the simulation. Keeps track of arguments passed,
//reading in data, and exporting data. Arguments passed is handled by specArgs class,
//reading in and constructing asteroids by inputFileHandler class, and writing to
//file is handled by outputFileHandler.

dataHandler::dataHandler(int argc, char** argv)
{
	allIsWell = true;
	argsHandler = specArgs(argc, argv);
	if (isGood())
		reader = inputFileHandler(argsHandler.getRandCount());
}

dataHandler::dataHandler(){}

void dataHandler::handOffResults(
	std::vector<vec3d> inputRelData,
	std::vector<std::vector<vec3d>> inputAbsData,
	std::vector<std::vector<vec3d>> inputUnitVectors,
	bool successfulSimulation){
	writer = outputFileHandler(inputRelData, inputAbsData, inputUnitVectors);
	allIsWell = successfulSimulation;
}

bool dataHandler::writeToFile(){
	#ifdef DEBUG_MODE
	print("writing to file");
	#endif
	writer.writeDataToFile(
		argsHandler.getFileName(
			reader.getSeed()),
		generateHeaderString(
			reader.getSeed(),
			argsHandler.getSelection(),
			simulationBodies,
			argsHandler.generateSampleTimes(),
			epochToYear(
				argsHandler.getInitialEpoch()),
			argsHandler.getComment()), 
		simTimeToFileTime(
			argsHandler.generateSampleTimes(),
			argsHandler.getInitialEpoch()),
		reader.getSeed());
}

vector<Orbit> dataHandler::generateAsteroids()
{
	#ifdef DEBUG_MODE
	print("generating asteroids");
	#endif
	vector<int> selection = argsHandler.getSelection();
	simulationBodies = vector<Orbit>();
	float initialEpoch = argsHandler.getInitialEpoch();
	for (int i = 0; i < selection.size(); i++){
		simulationBodies.push_back(reader.generateOrbit(selection[i], initialEpoch));
	}
	return simulationBodies;
}

vector<float> dataHandler::generateSampleTimes()
{
	return argsHandler.generateSampleTimes();
}

vector<double> dataHandler::generateFileTimes(){
	return simTimeToFileTime(generateSampleTimes(), argsHandler.getInitialEpoch());
}

LISA dataHandler::generateLISA()
{
	return LISA();
}

string dataHandler::generateFileName()
{
	return argsHandler.getFileName(reader.getSeed());
}

bool dataHandler::isGood()
{
	return argsHandler.isGood() && allIsWell;
}

void dataHandler::killProgram(){
	allIsWell = false;
}
