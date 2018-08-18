#include "globalSwitchStatements.h"
#include "inputFileHandler.h"

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>

#include "Orbit.cuh"

using namespace std;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "input: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define FLOATS_PER_ORBIT 20
#define SEMI_MAJOR 0
#define ECCENTRICITY 1
#define INCLINE 2
#define LONG_ASCEND_NODE 3
#define ARG_PERI 4
#define MEAN_ANOMALY 5
#define SIGMA_SEMI_MAJOR 6
#define SIGMA_ECC 7
#define SIGMA_INCLINE 8
#define SIGMA_LONG_ASCEND 9
#define SIGMA_ARG_PERI 10
#define SIGMA_MEAN_ANOM 11
#define GM 12
#define MG 12
#define ABSOLUTE_MAGNITUDE 13
#define SIGMA_ABS_MAG 14
#define ALBEDO 15
#define TYPE 16
#define TYPE_SIMPLE 17
#define SIMPLE_TYPE 17
#define EPOCH_NORMAL 18
#define EPOCH_MJD 19

#define BASE_FILE_LOC string("../../Reference_File/Asteroid_Data/")//R"(C:\Users\David\Documents\reference file\asteroid data\)"
#define NAMES string("Names.csv")
#define PARAMS string("Orbital_Parameters.csv")
#define MISC string("Misc_Specs.csv")
#define EPOCHS string("Epochs.csv")
#define RAW_BYTE_DATA string("Orbital_Data.bin")

#include "defineOrbitalValues.h"

#define NORM_CONST 0.39894228f
#define METRO_ALGORITHM_DX 3

//reads in and generates asteroids. Initializes by reading asteroid data and generating
//gaussian distribution. Can then generate vector<Orbit> of asteroids when given
//a vector<int> of indices to generate. findMass has important information.

float random(float min, float max)
{
	int r = rand();
	float firstRand = (float)r / RAND_MAX;
	return firstRand * (max - min) + min;
}

float random(float max)
{
	return random(0, max);
}

vector<string> readCSVtoStringVector(string location)
{
	ifstream file(location);
	vector<string> output;
	string contents{ istreambuf_iterator<char>(file), istreambuf_iterator<char>() };
	vector<int> breakPoints;
	int i = 0;
	for (string::iterator it = contents.begin(); it != contents.end(); ++it)
	{
		if (*it == ',' || *it == '\n')
		{
			breakPoints.push_back(i);
		}
		i++;
	}
	for (int i = 1; i < breakPoints.size(); i++){
		string temp = contents.substr(breakPoints[i - 1] + 1,
			(breakPoints[i] - breakPoints[i - 1] - 2));
		output.push_back(temp);
	}
	return output;
}

float normalDistribution(float x)
{
	return NORM_CONST * expf(-x*x / 2);
}

vector<float> metroAlg(float(*distributionFunction)(float),
	int outputCount, float initialPoint, int dryRunCount, float dx)
{
	vector<float> output;
	float x = initialPoint;
	float y = distributionFunction(x);
	for (int i = 0; i < outputCount + dryRunCount; i++)
	{
		float nextX = x + dx * ((float)rand() * 2 / RAND_MAX - 1);
		float nextY = distributionFunction(nextX);
		float weight = nextY / y;
		if (weight >= 1)
		{
			x = nextX;
			y = nextY;
		}
		else
		{
			if (weight * RAND_MAX > rand())
			{
				x = nextX;
				y = nextY;
			}
		}
		if (i > dryRunCount)
		{
			output.push_back(x);
		}
	}
	return output;
}

vector<float> readRawFloatFile(string location)
{
	ifstream fin(location, ios::binary);
	vector<char> buffer = vector<char>(istreambuf_iterator<char>(fin), istreambuf_iterator<char>());
	int length = buffer.size() / sizeof(float);
	float* temp = (float*)buffer.data();
	vector<float> output{ temp, temp + length };
	return output;
}

float randType()
{
	int r = rand();
	if (r < .75 * RAND_MAX) return 0;
	if (r < .92 * RAND_MAX) return 1;
	else return 2;
}



inputFileHandler::inputFileHandler()
{
}

inputFileHandler::inputFileHandler(int numberOfRandomNumbers)
{
	seed = (int)time(0);
	srand(seed);
	orbitalData = readRawFloatFile(BASE_FILE_LOC + RAW_BYTE_DATA);
	//dataLength = orbitalData.size() / FLOATS_PER_ORBIT;
	normalRands = metroAlg(normalDistribution, numberOfRandomNumbers, 0, 1000, METRO_ALGORITHM_DX);
	currentRandIndex = 0;
	orbitNames = readCSVtoStringVector(BASE_FILE_LOC + NAMES);
}

vector<float> inputFileHandler::getOrbitData(int index)
{
	vector<float> output;
	for (int i = FLOATS_PER_ORBIT * index; i < FLOATS_PER_ORBIT * (index + 1); i++)
	{
		output.push_back(orbitalData[i]);
	}
	return output;
}

float inputFileHandler::findMass(float absMag, float type, float albedo)
{
	//when mass is not given, this method will generate an estimate
	//based on albedo, type, and absolute magnitude.
	if (type == -1)//type not specified
	{
		if (albedo == 0)//albedo not specified
		{
			float x, dx;
			type = randType();//random type based on current predicted compositions
			switch ((int)type)
			{
			case (0)://carbonaceous
				x = -2.9046f;
				dx = 0.602f;
				break;
			case(1)://stony
				x = -1.908;
				dx = 0.394;
				break;
			case(2)://metalic
				x = -1.956;
				dx = 0.3466;
				break;
			}
			albedo = expf(getRandValue(x, dx));
		}
		else//no type but albedo given. guess type based on albedo
		{
			if (albedo < .8) type = 0;
			else if (albedo < 1.3) type = randType();
			else if (albedo < 2.3)
			{
				if (rand() < .68 * RAND_MAX) type = 1;
				else type = 2;
			}
			else type = 1;
		}
	}
	else if (albedo == 0)//type specified but albedo not specified
	{//same code as in first case, just without randomized type
		float x, dx;
		switch ((int)type)
		{
		case (0):
			x = -2.9046f;
			dx = 0.602f;
			break;
		case(1):
			x = -1.908;
			dx = 0.394;
			break;
		case(2):
			x = -1.956;
			dx = 0.3466;
			break;
		}
		albedo = expf(getRandValue(x, dx));
	}//if both specified then we don't need to estimate anything above
	float density;
	switch ((int)type)
	{
	case(0)://carbonaceous density
		density = 1380;
		break;
	case(1)://stony density
		density = 2710;
		break;
	case(2)://metalic density
		density = 5320;
		break;
	}
	//take log to linear for absMag in middle statement
	//and handle the fact that we are working with area and so there is a sqrt
	//in there to get to radius. specifics of numbers were taken from wikipedia.
	float radius = 1329000 * powf(10, -absMag / 5) / (2 * sqrtf(albedo));
	//basic volume calculation
	float volume = 4 * PI / 3 * radius * radius * radius;
	return volume * density;//basic mass calculation
}

float inputFileHandler::getRandValue(float mean, float standardDeviation)
{
	//all number were generated in monty-carlo loop at object instantiation.
	//now we just sample the generated distribution. This is done by taking
	//a large step size through the list. The step size is prime to avoid
	//having it be a factor of the total size (since this would result in
	//missing most of the values). Once done, the standard deviation we want
	//is achieved by multiplying by it and the mean by adding it.
	float x = normalRands[currentRandIndex];
	currentRandIndex += 127;
	if (currentRandIndex >= normalRands.size())
		currentRandIndex -= normalRands.size();
	x *= standardDeviation;
	x += mean;
	return x;
}

Orbit inputFileHandler::generateOrbit(int index, float time0)
{
	//The macros at the top specify what position in the input file each datum of
	//a given orbit is in. So here, we can specify the position by name.
	vector<float> orbitData = getOrbitData(index);
	float absMag;
	if (orbitData[SIGMA_ABS_MAG] == 0)//if no error bars for absolute magnitude are given,
		//take it to be .2.
		absMag = getRandValue(orbitData[ABSOLUTE_MAGNITUDE], .2f);
	else absMag = getRandValue(orbitData[ABSOLUTE_MAGNITUDE], orbitData[SIGMA_ABS_MAG]);
	float mass;
	//if mass is not given, the use findMass function to create an estimate.
	//otherwise, use the given mass (with unit conversions to get into SI).
	if (orbitData[MG] == 0) mass = findMass(absMag, orbitData[SIMPLE_TYPE], orbitData[ALBEDO]);
	else mass = orbitData[MG] / G * 1e9f;
	return Orbit(getRandValue(orbitData[MEAN_ANOMALY], orbitData[SIGMA_MEAN_ANOM]),
		getRandValue(orbitData[ECCENTRICITY], orbitData[SIGMA_ECC]),
		AU*getRandValue(orbitData[SEMI_MAJOR], orbitData[SIGMA_SEMI_MAJOR]),
		getRandValue(orbitData[INCLINE], orbitData[SIGMA_INCLINE]),
		getRandValue(orbitData[LONG_ASCEND_NODE], orbitData[SIGMA_LONG_ASCEND]),
		getRandValue(orbitData[ARG_PERI], orbitData[SIGMA_ARG_PERI]),
		false, orbitalData[EPOCH_NORMAL], time0, mass, getName(index));
}

vector<Orbit> inputFileHandler::generateOrbits(vector<int> indices, float time0)
{
	vector<Orbit> output;
	for (vector<int>::iterator it = indices.begin(); it != indices.end(); ++it)
	{
		output.push_back(generateOrbit(*it, time0));
	}
	return output;
}

string inputFileHandler::getName(int index)
{
	return orbitNames[index];
}

unsigned int inputFileHandler::getSeed()
{
	return seed;
}