#include "globalSwitchStatements.h"
//#include <stdio.h>
#include <iostream>
#include <math.h>
#include <chrono>
#include <string>
#include <fstream>
#include <thread>

//#include "Orbit.cuh"
//#include "vec3d.cuh"
#include "kernel.cuh"
//#include "fileHandler.h"
#include "dataHandler.h"
// #include "inputDataHandler.h"

using namespace std;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "sim: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define NORMAL_RANDOM_NUMBER_COUNT 100000

#define HALFPI 1.5707963267949f
#define TWOPI 6.2831853071796f
#define PI 3.1415926535898f
#define au 149597870700.0f
#define year 31536000
#define armLength 5.0e9f
#define G 6.67408e-11f

int dynamicExecution(int argc, char** argv)
{
	print("Initializing.");
	dataHandler data = dataHandler(argc, argv);
	if (data.isGood())
	{
		print("Initialization Complete. Beginning Simulation");
		clock_t time0 = clock();
		bool isGood = executeCudaCode(data);
		clock_t time1 = clock();
		float ellapsedTime = (float)(time1 - time0) / CLOCKS_PER_SEC;
		print("Simulation Completed in " << ellapsedTime << " seconds.");
		if (isGood)
		{
			data.writeToFile();
		}
		else
		{
			print("Critical Failure in Simulation!");
		}
	}
	else
	{
		print("Bad Initialization.");
	}

	return 0;
}

int main(int argc, char** argv)
{
	return dynamicExecution(argc, argv);
	//return debugExecution();
}
