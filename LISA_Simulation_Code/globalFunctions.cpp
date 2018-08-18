#include "globalSwitchStatements.h"
#include "globalFunctions.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>

#include "Orbit.cuh"
#include "vec3d.cuh"

using namespace std;
// using namespace arguments;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "global: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#include "defineOrbitalValues.h"

float lisa::yearToEpoch(float inputYear){
	inputYear -= 2030;
	inputYear *= 365;
	inputYear += 62502;
	return inputYear;
}
float lisa::epochToYear(float inputEpoch){
	inputEpoch -= 62502;
	inputEpoch /= 365;
	inputEpoch += 2030;
	return inputEpoch;
}

vector<double> lisa::simTimeToFileTime(vector<float> simTime, float initialEpoch){
	#ifdef DEBUG_MODE
	print("translating sim time to file time");
	#endif
	vector<double> output;
	double initialYear = epochToYear(initialEpoch);
	for (int i = 0; i < simTime.size(); i++){
		output.push_back((double)simTime[i] / YEAR + initialYear);
	}
	return output;
}

string lisa::generateHeaderString(
	unsigned int seed,
	vector<int> asteroidIndices,
	vector<Orbit> bodiesSimulated,
	vector<float> timesSampled,
	float initialYear,
	string comment)
{
	#ifdef DEBUG_MODE
	print("generating file header");
	#endif
	string output = "date=";
	//generate todays date
	time_t t = time(0);
    tm now = *(localtime(&t));
    int year = now.tm_year + 1900;
    int month = now.tm_mon + 1;
    int day = now.tm_mday;
    string date = to_string(year);
    if (month < 10){
    	date += "0";
    }
    date += to_string(month);
    if (day < 10){
    	date += "0";
    }
    date += to_string(day);
    //seed
    output += date + "\nseed=";
    output += to_string(seed) + "\ninitial year=";
    //initial year and generate final year
    output += to_string(initialYear) + "\nfinal year=";
    float dt = timesSampled[1] - timesSampled[0];
    float lastTime = timesSampled.back() + dt;
    float endYear = initialYear + lastTime / year;
    output += to_string(endYear) + "\ntime between samples=";
    output += to_string(dt) + "\nnumber of time samples=";
    output += to_string(timesSampled.size()) + "\nfirst asteroid index=";
    output += to_string(asteroidIndices[0]) + "\nlast asteroid index=";
    output += to_string(asteroidIndices.back()) + "\ntotal number of asteroids=";
    output += to_string(bodiesSimulated.size()) + "\nfirst 10 asteroid names=";
    int temp = 10;
    output += bodiesSimulated[0].getName();
    if (bodiesSimulated.size() < 10) temp = bodiesSimulated.size();
    for (int i = 1; i < temp; i++){
    	output += "," + bodiesSimulated[i].getName();
    }
    output += "\ncomments=";
    output += comment;
	return output;
}

vector<double> lisa::floatVectorToDoubleVector(vector<float> input){
    vector<double> output = vector<double>();
    for (int i = 0; i < input.size(); i++){
        output.push_back(input[i]);
    }
    return output;
}

vector<double> lisa::vecVectorToDoubleVector(vector<cudaUtil::vec3d> input){
    vector<double> output = vector<double>();
    for (int i = 0; i < input.size(); i++){
        output = concat(output, floatVectorToDoubleVector(input[i].toFloatVector()));
    }
    return output;
}

void writeCSVfile(vector<vector<double>> data, string location){
    #ifdef DEBUG_MODE
    print("writing csv file")
    #endif
    ofstream fout;
    fout.open(location, ios::out | ios::trunc);
    fout << to_string(data[0][0]);
    for (int j = 1; j < data[0].size(); j++)
    {
        fout << "," << data[0][j];
    }
    for (int i = 1; i < data.size(); i++)
    {
        fout << "\n";
        fout << to_string(data[i][0]);
        for (int j = 1; j < data[i].size(); j++)
        {
            fout << "," << data[i][j];
        }
    }
    fout.close();
}

string str(double input){
    double exponent = log10(input);
    int intExp;
    if (exponent >= 0){
        intExp = int(exponent);
    }
    else{
        intExp = int(exponent) - 1;
    }
    
    string output = to_string(input / pow(10, intExp));
    if (intExp != 0){
        output += "e";
        output += to_string(intExp);
    }
    return output;
}