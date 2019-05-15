#include "globalSwitchStatements.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>

#include "dataHandler.h"
#include "vec3d.cuh"
#include "Orbit.cuh"
#include "globalFunctions.h"

using namespace std;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "mass: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define BASE_FILE_LOC string("../../Reference_File/Data_Dumps/")
#define MASS_ANALYSIS string("Mass_Analysis/")
#define SUFFIX ".csv"

#define foreach(type, name, inputVector)\
	for (vector<type>::iterator name = inputVector.begin();\
	name != inputVector.end(); name++)

///////////////NASA known masses////////////////////////////////////////////////
		
// vector<int> interestingAsteroidIndices = vector<int>({1, 2, 3, 4, 21,
// 	45, 140, 216, 243, 253, 433, 1566, 1620, 1862, 2060, 2530, 2703,
// 	2867, 3352, 3840, 4179, 4769, 4979, 5535, 9969, 25143, 101955});

// vector<float> acceptedMasses = vector<float>({939300, 205000, 20000,
// 	259000, 1700, 6100, 1500, 0, 100, 103.3, 6.69, .001, .004, .002, 4000,
// 	0, 0, 0, 0, 0, .05, .0005, .2, 0, 0, .000035, .00014});

//////////////////////////////////////////////////////////////////////////////

/////////////////////////Wikipedia!!///////////////////////////////////////////

// vector<int> interestingAsteroidIndices = vector<int>({1, 4, 2, 10, 31,
// 	704, 511, 532, 15, 3, 16, 52, 88, 7, 13, 423, 29, 87, 48});

// vector<float> acceptedMasses = vector<float>({939300, 259076, 201000,
// 	86700, 58100, 38800, 37700, 33000, 31800, 28600, 22700, 18300,
// 	16200, 16000, 16000, 15200, 14780, 12000});

///////////////////////////////////////////////////////////////////////////////////

//////////////////////////EAR Compilation////////////////////////////////////////

vector<int> interestingAsteroidIndices = vector<int>({1, 2, 3, 4, 6, 7, 8,
	9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 29, 31, 41, 45, 47, 48,
	49, 52, 65, 87, 88, 90, 107, 121, 130, 189, 243, 253, 283, 324, 379, 433,
	444, 451, 511, 532, 617, 702, 704, 762, 804, 3749, 25143, 26308, 42355,
	47171, 50000, 55636, 58534, 65489, 66391, 66652, 88611, 90482, 136199, 136108,
	134860, 296962, 216});

vector<float> acceptedMasses = vector<float>({});

/////////////////////////////////////////////////////////////////////////////////

int main()
{
	//remove zero mass values (for first dataset)
	for (int i = acceptedMasses.size() - 1; i >= 0; i--){
		if (acceptedMasses[i] == 0){
			acceptedMasses.erase(acceptedMasses.begin() + i);
			interestingAsteroidIndices.erase(interestingAsteroidIndices.begin() + i);
		}
	}

	//produce arguments to feed into dataHandler constructor
	char** inducingArgs = new char*[acceptedMasses.size() + 2];
	inducingArgs[0] = new char[string("mass.exe").length()];
	strcpy(inducingArgs[0], string("mass.exe").c_str());
	inducingArgs[1] = new char[string("-s").length()];
	strcpy(inducingArgs[1], string("-s").c_str());
	for (int i = 0; i < acceptedMasses.size(); i++){
		inducingArgs[i + 2] = new char[to_string(interestingAsteroidIndices[i] - 1).length()];
		strcpy(inducingArgs[i + 2], to_string(interestingAsteroidIndices[i] - 1).c_str());
	}
	
	dataHandler handler = dataHandler(acceptedMasses.size() + 2, inducingArgs);

	//grab asteroids
	vector<Orbit> asteroids = handler.generateAsteroids();
	vector<double> masses = vector<double>();
	
	//grab masses from asteroids
	for (int i = 0; i < asteroids.size(); i++){
		masses.push_back(asteroids[i].getMass());
	}
	
	//construct list that writeCSVfile can handle
	vector<vector<double>> writingList = vector<vector<double>>();
	
	for (int i = 0; i < masses.size(); i++){
		vector<double> tempList = vector<double>();
		tempList.push_back(interestingAsteroidIndices[i]);
		tempList.push_back(acceptedMasses[i]);
		tempList.push_back(masses[i] / pow(10, 15));
		writingList.push_back(tempList);
	}
	
	writeCSVfile(writingList, BASE_FILE_LOC + MASS_ANALYSIS + string("spec") + SUFFIX);

	return 0;
}
