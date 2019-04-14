#include "globalSwitchStatements.h"
#include "dataHandler.h"

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>

#include "Orbit.h"
#include "LISA.h"
#include "specArgs.h"
#include "inputFileHandler.h"
#include "globalFunctions.h"

#include "defineOrbitalValues.h"

using namespace std;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "data: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

typedef vsu::ProductSpace<double, double, double> vec3;

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

LISA dataHandler::generateLISA()
{
	const double secondsPerDay = 24 * 3600;
	// const double secondsPerDay = 1;
	const double orbitalParameterEpoch = 28324.75;
	// const double epoch2030 = 62502;
	// const double epochOffset = -365;
	// const double startingEpoch = epoch2030 + epochOffset;
	double startingEpoch = argsHandler.getInitialEpoch();
	print("initial epoch: " << startingEpoch);
	// earth = LISAutils.orbit(math.radians(200.7), 0.01671, math.radians(0),
	// 	math.radians(-11.261), math.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
	// 	simTime = startEpoch, paramTime = 58324.75*24*3600, name = 'earth')
	Orbit earth(radians(200.7),//mean anomaly
		0.01671,//ecc
		AU,//semi major
		radians(0),//incline
		radians(-11.261),//longitude ascending
		radians(114.2078),//argument perihelion
		true,//provided mean anomaly, not true anomaly
		orbitalParameterEpoch * secondsPerDay,//parameter epoch
		startingEpoch * secondsPerDay,//initial simulation epoch
		5.972e24,//Earth mass
		"Earth");//name
	return LISA(-PI * 2 / 3, atan2(earth.y(), earth.x()) - PI * 20 / 180);
}

bool dataHandler::isGood()
{
	return argsHandler.isGood() && allIsWell;
}
