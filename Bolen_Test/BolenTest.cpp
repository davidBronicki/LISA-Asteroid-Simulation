
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "inputFileHandler.h"
#include "globalFunctions.h"

#include "defineOrbitalValues.h"

using namespace std;
using namespace vsu;
using namespace lisa;



int main()
{
	inputFileHandler dataHandler(1000000);// initialize data

	vector<int> indices;
	for (int i = 0; i < 10; ++i)//first 10 asteroids
	{
		indices.push_back(i);
	}

	double startingEpoch = yearToEpoch(2025);

	vector<Orbit> asteroids = dataHandler.generateOrbits(indices, startingEpoch * 24 * 3600);

	vector<vector<double>> outputData;//list of each asteroid's x, y, z, and mass for each point of time

	for (double t = 0; t < 5 * 365 * 24 * 3600; t += 24 * 3600)
	{
		//running for 5 years, and dt = 1 day.
		vector<double> dataForThisTimeStep;
		//dataForThisTimeStep = []

		for (auto&& asteroid : asteroids)
		{
			asteroid.setTime(t);

			dataForThisTimeStep.push_back(asteroid.x());
			dataForThisTimeStep.push_back(asteroid.y());
			dataForThisTimeStep.push_back(asteroid.z());
			dataForThisTimeStep.push_back(asteroid.getMass());
		}

		outputData.push_back(dataForThisTimeStep);
		//outputData.append(dataForThisTimeStep)
	}


	vector<string> columnHeaders;

	for (auto&& asteroid : asteroids)
	{
		columnHeaders.push_back(asteroid.getName() + "__x");
		columnHeaders.push_back(asteroid.getName() + "__y");
		columnHeaders.push_back(asteroid.getName() + "__z");
		columnHeaders.push_back(asteroid.getName() + "__mass");
	}


	writeCSVfile(outputData,
		columnHeaders,
		"Bolen_Test_Output_Directory/Bolen_Sample_Data.csv");
}

