#include "globalSwitchStatements.h"
#include "specArgs.h"

#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <math.h>

#include "defineOrbitalValues.h"

using namespace std;
using namespace lisa;
using namespace arguments;

#ifdef DEBUG_MODE
#define print(input) cout << "spec: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define DEFAULT_RAND_COUNT 100000
#define DEFAULT_TIME_COUNT 10000
#define DEFAULT_AST_COUNT 10000

//arguments handler specific to this simulation. sets up a number of available arguments
//to pass (see help string). Uses argsHandler class to fascilitate. Once instantiated,
//it can construct key data such as sample times, initial epoch, and the indices of
//asteroids selected.

specArgs::specArgs(int argc, char** argv)
{
	// killProgram = false;
	handler = argumentsHandler();
	handler.addFlag(vector<string>({"n", "name"}), STRING, -1);
	handler.addFlag(vector<string>({"s", "select"}), INT, -1);
	handler.addFlag(vector<string>({"r", "range"}), INT, 2);
	handler.addFlag(vector<string>({"c", "count"}), INT, 1);
	handler.addFlag(vector<string>({"rc", "randcount"}), INT, 1);
	handler.addFlag(vector<string>({"tc", "timecount"}), INT, 1);
	handler.addFlag(vector<string>({"st", "ts", "sampletime", "samplerate"}), FLOAT, 1);
	handler.addFlag(vector<string>({"tr", "timerange"}), FLOAT, 2);
	handler.addFlag(vector<string>({"comment", "comments"}), STRING, -1);
	handler.setHelp(
R"(-n/-name to name the file
-s/-select with series of integers to select those asteroids
-r/-range with two integers to select range (excluding second number)
-c/-count with one integer to select first n asteroids
-rc/-randcount with one integer to select number of random numbers to generate
-tc/-timecount with one integer to select number of sample points in time
-sr/-samplerate with one float to select sample rate
-tr/-timerange with two float to select start and stop times in years
-comment with group of strings to say what you think about the simulation)");
	handler.passArgs(argc, argv);
}

specArgs::specArgs(){}

string specArgs::getComment(){
	#ifdef DEBUG_MODE
	print("getting comment");
	#endif
	if (handler.isFlagCalled("comment")){
		vector<string> temp = *(vector<string>*)(handler.getFlagValue("comment"));
		string output = "";
		output += temp[0];
		for (int i = 1; i < temp.size(); i++){
			output += " " + temp[i];
		}
		return output;
	}
	else{
		return "N/A";
	}
}

vector<int> specArgs::getSelection(){
	#ifdef DEBUG_MODE
	print("getting selection indices");
	#endif
	vector<int> selection = vector<int>();
	if (handler.isFlagCalled("count")){
		vector<int> temp = *(vector<int>*)(handler.getFlagValue("count"));
		for (int i = 0; i < temp[0]; i++){
			selection.push_back(i);
		}
	}
	else if (handler.isFlagCalled("range")){
		vector<int> temp = *(vector<int>*)(handler.getFlagValue("range"));
		for (int i = temp[0]; i < temp[1]; i++){
			selection.push_back(i);
		}
	}
	else if (handler.isFlagCalled("select")){
		selection = *(vector<int>*)(handler.getFlagValue("select"));
	}
	else{
		for (int i = 0; i < DEFAULT_AST_COUNT; i++){
			selection.push_back(i);
		}
	}
	return selection;
}

vector<float> specArgs::generateSampleTimes(){
	#ifdef DEBUG_MODE
	print("generating sim times");
	#endif
	vector<float> output;
	float t0;
	float t1;
	if (handler.isFlagCalled("timerange")){
		vector<double> temp = *(vector<double>*)(handler.getFlagValue("timerange"));
		t0 = min(temp[0], temp[1]);
		t1 = max(temp[0], temp[1]);
	}
	else{
		t0 = 0;
		t1 = 1;
	}
	t0 *= YEAR;
	t1 *= YEAR;
	float dt;
	int sampleTimeCount;
	if (handler.isFlagCalled("timecount")){
		vector<int> temp = *(vector<int>*)(handler.getFlagValue("timecount"));
		dt = (t1 - t0) / temp[0];
		sampleTimeCount = temp[0];
	}
	else if (handler.isFlagCalled("samplerate")){
		vector<float> temp = *(vector<float>*)(handler.getFlagValue("samplerate"));
		dt = temp[0];
		sampleTimeCount = (t1 - t0) / dt;
	}
	else{
		dt = (t1 - t0) / DEFAULT_TIME_COUNT;
		sampleTimeCount = DEFAULT_TIME_COUNT;
	}
	for (int i = 0; i < sampleTimeCount; i++){
		output.push_back(i * dt);
	}
	return output;
}

float specArgs::getInitialEpoch(){
	#ifdef DEBUG_MODE
	print("getting initial epoch");
	#endif
	if (handler.isFlagCalled("timerange")){
		vector<int> temp = *(vector<int>*)(handler.getFlagValue("timerange"));
		return min(temp[0], temp[1]);
	}
	else{
		return 0;
	}
}

string specArgs::getFileName(unsigned int seed){
	#ifdef DEBUG_MODE
	print("getting file name");
	#endif
	if (handler.isFlagCalled("name")){
		string output = "";
		vector<string> temp = *(vector<string>*)(handler.getFlagValue("name"));
		output += temp[0];
		for (int i = 1; i < temp.size(); i++){
			output += " " + temp[i];
		}
		return output;
	}
	else{
		print("File name not provided. Defaulting to seed.");
		return to_string(seed);
	}
}

int specArgs::getRandCount(){
	#ifdef DEBUG_MODE
	print("getting rand count");
	#endif
	if (handler.isFlagCalled("randcount")){
		vector<int> temp = *(vector<int>*)(handler.getFlagValue("randcount"));
		return temp[0];
	}
	else{
		return DEFAULT_RAND_COUNT;
	}
}

bool specArgs::isGood(){
	#ifdef DEBUG_MODE
	print("checking status");
	#endif
	int temp = 0;
	temp += handler.isFlagCalled("count");
	temp += handler.isFlagCalled("range");
	temp += handler.isFlagCalled("select");
	if (temp > 1) return false;
	temp = 0;
	temp += handler.isFlagCalled("timecount");
	temp += handler.isFlagCalled("samplerate");
	if (temp > 1) return false;
	return !handler.isHelpCalled();
}
