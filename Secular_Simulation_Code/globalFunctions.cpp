#include "globalSwitchStatements.h"
#include "globalFunctions.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>

#include "Orbit.h"
#include "VectorSpace.h"

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

// vector<double> lisa::vecVectorToDoubleVector(vector<cudaUtil::vec3d> input){
//     vector<double> output = vector<double>();
//     for (int i = 0; i < input.size(); i++){
//         output = concat(output, floatVectorToDoubleVector(input[i].toFloatVector()));
//     }
//     return output;
// }

void writeCSVfile(vector<vector<double>> data, string location){
    #ifdef DEBUG_MODE
    print("writing csv file");
    #endif
    ofstream fout;
    fout.open(location, ios::out | ios::trunc);
    fout.precision(16);
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

#include "undefineOrbitalValues.h"
