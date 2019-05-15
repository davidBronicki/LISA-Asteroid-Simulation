#include "globalSwitchStatements.h"
#include "outputFileHandler.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "vec3d.cuh"
#include "globalFunctions.h"
//#include "fileHandler.h"

using namespace std;
using namespace lisa;
using namespace cudaUtil;

#ifdef DEBUG_MODE
#define print(input) cout << "output: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define BASE_FILE_LOC string("../../Reference_File/Data_Dumps/")
#define RELATIVE string("Arm_Acceleration/")
#define ABSOLUTE_ECLIPTIC string("Ecliptic_Acceleration/")
#define UNIT_VECTOR string("LISA_Unit_Vector/")
#define HEADER string("Header/")
#define SUFFIX ".csv"

outputFileHandler::outputFileHandler(){}

outputFileHandler::outputFileHandler(
	vector<vec3d> inputRelData,
	vector<vector<vec3d>> inputAbsData,
	vector<vector<vec3d>> inputUnitVectors)
{
	relativeAccelData = inputRelData;
	absoluteAccelData = inputAbsData;
	unitVectorData = inputUnitVectors;
}

void writeHeader(string headerString, string location){
	ofstream fout;
	fout.open(location, ios::out | ios::trunc);
	fout << headerString;
	fout.close();
}

void outputFileHandler::writeDataToFile(
	string fileName,
	string headerString,
	vector<double> fileTimes,
	unsigned int seed)
{
	#ifdef DEBUG_MODE
	print("writing to file");
	#endif
	vector<vector<double>> relativeAccelWritingData = vector<vector<double>>();
	vector<vector<double>> absoluteEclipticAccelWritingData = vector<vector<double>>();
	vector<vector<double>> unitVectorWritingData = vector<vector<double>>();
	for (int i = 0; i < fileTimes.size(); i++)
	{
		vector<double> tempRelative = vector<double>({fileTimes[i]});
		tempRelative = concat(tempRelative,
			floatVectorToDoubleVector(relativeAccelData[i].toFloatVector()));
		relativeAccelWritingData.push_back(tempRelative);

		vector<double> tempEcliptic = vector<double>({fileTimes[i]});
		tempEcliptic = concat(tempEcliptic, vecVectorToDoubleVector(absoluteAccelData[i]));
		absoluteEclipticAccelWritingData.push_back(tempEcliptic);

		vector<double> tempUnit = vector<double>({fileTimes[i]});
		tempUnit = concat(tempUnit, vecVectorToDoubleVector(unitVectorData[i]));
		unitVectorWritingData.push_back(tempUnit);
	}
	writeCSVfile(relativeAccelWritingData, BASE_FILE_LOC + RELATIVE + fileName + SUFFIX);
	writeCSVfile(absoluteEclipticAccelWritingData, BASE_FILE_LOC + ABSOLUTE_ECLIPTIC + fileName + SUFFIX);
	writeCSVfile(unitVectorWritingData, BASE_FILE_LOC + UNIT_VECTOR + fileName + SUFFIX);
	writeHeader(headerString, BASE_FILE_LOC + HEADER + fileName);
}
