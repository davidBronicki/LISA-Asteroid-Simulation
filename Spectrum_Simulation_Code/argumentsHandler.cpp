#include "globalSwitchStatements.h"
#include "argumentsHandler.h"

#include <tuple>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace arguments;

#ifdef DEBUG_MODE
#define print(input) cout << "args: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define DATA_TYPE int

#define INT 0
#define FLOAT 1
#define STRING 2
#define BOOL 3

//generic arguments handler. flags are added by specifying a list of names
//(vector<string>), data type (see #define's above), and data size (int).
//data size is how many arguments should be passed to the given flag, -1
//means any number is valid. Also setHelp(string) will setup a help
//display if the user runs -help or -h. passArgs is where the args
//go. It takes the exact set sent to main. isFlagCalled returns bool
//for if the given flag was called in args. getFlagValue returns
//a void* pointer. The programmer must cast it to vector<data type>.
//if no argument was passed for the flag, getFlagValue will return NULL.

vector<tuple<int, string>> locateFlags(vector<string> args){
	vector<tuple<int, string>> output = vector<tuple<int, string>>();
	for (int i = 0; i < args.size(); i++){
		if (args[i][0] == '-'){
			string temp = args[i];
			temp.erase(0, 1);
			output.push_back(make_tuple(i, temp));
		}
	}
	return output;
}

argumentsHandler::argumentsHandler(){
	dataCalled = vector<bool>();
	data = vector<void*>();
	dataTypes = vector<int>();
	dataSizes = vector<int>();
	flags = vector<vector<string>>();
	helpCalled = false;
}

argumentsHandler& argumentsHandler::addFlag(
		std::vector<std::string> flagNames, int dataType, int dataSize){
	flags.push_back(flagNames);
	dataTypes.push_back(dataType);
	dataSizes.push_back(dataSize);
	data.push_back(NULL);
	dataCalled.push_back(false);
	return *this;
}

int argumentsHandler::getFlagIndex(string flag){
	for (int i = 0; i < flags.size(); i++){
		for (int j = 0; j < flags[i].size(); j++){
			if (flag == flags[i][j]) return i;
		}
	}
	print("flag not found. Flag = " << flag);
	return -1;
}

bool argumentsHandler::setFlagValue(string flag, vector<string> rawDatum){
	int index = getFlagIndex(flag);
	if (index == -1) {
		return false;
	}
	DATA_TYPE dType = dataTypes[index];
	int size = dataSizes[index];
	if (size >= 0 && size != rawDatum.size()){
		print("invalid number of arguments. Flag = " << flag);
		return false;
	}
	if (dType == INT){
		vector<int>* temp = new vector<int>();
		for (int i = 0; i < rawDatum.size(); i++){
			try{
				temp->push_back(stoi(rawDatum[i]));
			}
			catch(...){
				print("invalid int passed. Flag = " << flag);
				return false;
			}
		}
		data[index] = (void*)temp;
	}
	if (dType == FLOAT){
		vector<double>* temp = new vector<double>();
		for (int i = 0; i < rawDatum.size(); i++){
			try{
				temp->push_back(stod(rawDatum[i]));
			}
			catch(...){
				print("invalid float passed. Flag = " << flag);
				return false;
			}
		}
		data[index] = (void*)temp;
	}
	if (dType == STRING){
		vector<string>* temp = new vector<string>();
		for (int i = 0; i < rawDatum.size(); i++){
			temp->push_back(rawDatum[i]);
		}
		data[index] = (void*)temp;
	}
	if (dType == BOOL){
		vector<bool>* temp = new vector<bool>();
		for (int i = 0; i < rawDatum.size(); i++){
			if (rawDatum[i] == "true"
				|| rawDatum[i] == "t"
				|| rawDatum[i] == "1"
				|| rawDatum[i] == "True"
				|| rawDatum[i] == "T"
				|| rawDatum[i] == "TRUE")
				temp->push_back(true);
			else if (rawDatum[i] == "false"
				|| rawDatum[i] == "f"
				|| rawDatum[i] == "0"
				|| rawDatum[i] == "False"
				|| rawDatum[i] == "F"
				|| rawDatum[i] == "FALSE")
				temp->push_back(false);
			else{
				print("invalid boolean passed. Flag = " << flag);
				return false;
			}
		}
		data[index] = (void*)temp;
	}
	dataCalled[index] = true;
	return true;
}

void argumentsHandler::passArgs(int argc, char** argv){
	vector<string> args = vector<string>();
	for (int i = 1; i < argc; i++){
		args.push_back(string(argv[i]));
	}
	vector<tuple<int, string>> flagLocations = locateFlags(args);
	for (int i = 0; i < flagLocations.size(); i++){
		int index = get<0>(flagLocations[i]);
		string flag = get<1>(flagLocations[i]);
		if (flag == "h" || flag == "help"){
			print(helpString);
			helpCalled = true;
			continue;
		}
		int nextIndex;
		if (i != flagLocations.size() - 1){
			nextIndex = get<0>(flagLocations[i + 1]);
		}
		else{
			nextIndex = args.size();
		}
		vector<string> rawData;
		for (int j = index + 1; j < nextIndex; j++){
			rawData.push_back(args[j]);
		}
		setFlagValue(flag, rawData);
	}
}

void* argumentsHandler::getFlagValue(string flag){
	#ifdef DEBUG_MODE
	print("getting flag value");
	#endif
	int index = getFlagIndex(flag);
	if (index == -1) return NULL;
	return data[index];
}

bool argumentsHandler::isFlagCalled(string flag){
	int index = getFlagIndex(flag);
	if (index == -1) return false;
	return dataCalled[index];
}

argumentsHandler& argumentsHandler::setHelp(string inputHelpString){
	helpString = inputHelpString;
	return *this;
}

bool argumentsHandler::isHelpCalled(){
	return helpCalled;
}
