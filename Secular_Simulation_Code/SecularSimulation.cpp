#include "globalSwitchStatements.h"
#include <iostream>
#include <math.h>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>
#include <thread>

#include "LISA.h"
#include "globalFunctions.h"
#include "inputFileHandler.h"
#include "RungeKutta.h"

#include "defineOrbitalValues.h"

using namespace std;
using namespace vsu;
using namespace lisa;

//Some silly code for printing
//which can be disabled by macros
#ifdef DEBUG_MODE
template<typename... Ts>
void print(Ts... inputs){
	((cout << "sim: ") << ... << inputs) << endl;
}
// void print(){
// 	cout << endl;
// }
#else
#ifdef NO_PRINT
void print(T input){}
void print(){}
#else
void print(T input){
	cout << input << endl;
}
void print(){
	cout << endl;
}
#endif
#endif

/*
Define types for use in the simulation.
See VectorSpace.h for what this actually does
(or don't; its a very complicated piece of code.
These are essentially direct sums of vector spaces.)
*/
typedef ProductSpace<double, double, double> vec3;

typedef ProductSpace<vec3, vec3> state;
typedef ProductSpace<state, state, state> LISA_Vector;
typedef ProductSpace<double, LISA_Vector, LISA_Vector> totalState;

typedef ProductSpace<vec3, double> asteroidData;

vector<Orbit> asteroids;
vector<asteroidData> asteroidPositions;

//Inverse square law code
vec3 asteroidForce(const vec3& displacement, double mass){
	double distance = displacement.defaultMagnitude();
	return displacement * mass * G / pow(distance, 3);
}

//Specialized for the sun
vec3 centralForce(const vec3& position)
{
	return -1 * position * MU / pow(position.defaultMagnitude(), 3);
}

//apply sun force and asteroid forces to a single satellite
state oribitalDerivativeWithAsteroids(const state& input, double time){
	auto position = project<0>(input);
	auto velocity = project<1>(input);
	vec3 force;
	for (asteroidData asteroid : asteroidPositions)
	{
		force += asteroidForce(project<0>(asteroid) - position, project<1>(asteroid));
	}
	return state(velocity, force + centralForce(position));
}

//apply only sun force to a single satellite
state oribitalDerivativeWithoutAsteroids(const state& input, double time){
	auto position = project<0>(input);
	auto velocity = project<1>(input);
	return state(velocity, centralForce(position));
}

//apply forces to all LISA satellites
LISA_Vector LISA_DerivativeWithAsteroids(const LISA_Vector& input, double time){
	return LISA_Vector(oribitalDerivativeWithAsteroids(project<0>(input), time),
		oribitalDerivativeWithAsteroids(project<1>(input), time),
		oribitalDerivativeWithAsteroids(project<2>(input), time));
}

//apply sun force to all LISA satellites
LISA_Vector LISA_DerivativeWithoutAsteroids(const LISA_Vector& input, double time){
	return LISA_Vector(oribitalDerivativeWithoutAsteroids(project<0>(input), time),
		oribitalDerivativeWithoutAsteroids(project<1>(input), time),
		oribitalDerivativeWithoutAsteroids(project<2>(input), time));
}

//idk why this is defined here, but it's self explanatory
#define THREAD_COUNT 3

//apply forces to one LISA constellation with sun and asteroids (the perturbed LISA)
//and apply only the sun force to another constellation (the unperturbed LISA)
totalState dualDerivative(const totalState& input){
	double time = project<0>(input);
	vector<thread> threads;
	for (int i = 0; i < THREAD_COUNT; i++){
		threads.push_back(thread([](int threadNumber, double time){
				for (int j = threadNumber; j < asteroids.size(); j += THREAD_COUNT){
					asteroids[j].setTime(time);
					asteroidPositions[j] = asteroidData(asteroids[j].pos(), asteroids[j].getMass());
				}
			}, i , time));
	}
	for (int i = 0; i < THREAD_COUNT; i++){
		threads[i].join();
	}
	return totalState(1,
		LISA_DerivativeWithoutAsteroids(project<1>(input), time),
		LISA_DerivativeWithAsteroids(project<2>(input), time));
}

//These codes collect the error of a given time step. Used for
//the dynamic time step in Runge Kutta.
double stateSqrError(const vec3& currentState, const vec3& deltaState, double dt)
{
	return deltaState.defaultSquareMagnitude() / currentState.defaultSquareMagnitude();
}

double orbitalSqrError(const state& currentState, const state& deltaState, double dt)
{
	double e1 = stateSqrError(project<0>(currentState), project<0>(deltaState), dt);
	double e2 = stateSqrError(project<1>(currentState), project<1>(deltaState), dt);
	return e1 + e2;
}

double LISA_SqrError(const LISA_Vector& currentState, const LISA_Vector& deltaState, double dt)
{
	double e1 = orbitalSqrError(project<0>(currentState), project<0>(deltaState), dt);
	double e2 = orbitalSqrError(project<1>(currentState), project<1>(deltaState), dt);
	double e3 = orbitalSqrError(project<2>(currentState), project<2>(deltaState), dt);
	return e1 + e2 + e3;
}

double dualError(const totalState& currentState, const totalState& deltaState, double dt)
{
	double e1 = LISA_SqrError(project<1>(currentState), project<1>(deltaState), dt);
	double e2 = LISA_SqrError(project<2>(currentState), project<2>(deltaState), dt);
	return sqrt(e1 + e2);
}


//Converts an "orbit" object (with orbital elements)
//into a "state" object (with position and velocity)
state extractSat(Orbit input)
{
	return state(input.pos(), input.vel());
}

totalState extractLISA(LISA input)
{
	vector<Orbit> temp = input.getSats();
	LISA_Vector temp2(extractSat(temp[0]),
		extractSat(temp[1]),
		extractSat(temp[2]));
	return totalState(0, temp2, temp2);
}

//This is to help in the file writing procedure.
//Makes it so later we only need a single for loop.
vector<double> flattendLISA(LISA_Vector input)
{
	vector<double> output;
	vector<state> sats({project<0>(input),
		project<1>(input),
		project<2>(input)});
	vector<vec3> vecs;
	for (auto sat : sats)
	{
		vecs.push_back(project<0>(sat));
		vecs.push_back(project<1>(sat));
	}
	for (auto vec : vecs)
	{
		output.push_back(project<0>(vec));
		output.push_back(project<1>(vec));
		output.push_back(project<2>(vec));
	}
	return output;
}

//This is for printing the real time duration of the simulation
string formatTime(double time){
	int hrs = time / 3600;
	time -= hrs * 3600;
	int mins = time / 60;
	time -= mins * 60;
	int secs = time;
	string output;
	if (hrs != 0){
		output += to_string(hrs);
		output += ":";
	}
	if (hrs != 0 || mins != 0)
	{
		if (mins < 10){
			output += "0";
		}
		output += to_string(mins);
		output += ":";
	}
	if (secs < 10){
		output += "0";
	}
	output += to_string(secs);
	return output;
}

//info to put in header file (separate from the raw data file)
vector<string> generateHeaderStringList()
{
	vector<string> output;
	output.push_back("simulation_time");
	for (auto&& nameBase : vector<string>{"unperturbed_", "perturbed_"})
	{
		for (int satIndex = 0; satIndex < 3; ++satIndex)
		{
			for (auto&& vectorType : vector<string>{"_position_", "_velocity_"})
			{
				output.push_back(nameBase + to_string(satIndex) + vectorType + "x");
				output.push_back(nameBase + to_string(satIndex) + vectorType + "y");
				output.push_back(nameBase + to_string(satIndex) + vectorType + "z");
			}
		}
	}
	return output;
}

void runSimulation(const LISA& lisa, string simulationDescriptor)
{
	auto headerStringList = generateHeaderStringList();

	RungeKuttaSimulation<totalState> simulation(extractLISA(lisa), 1.0e4, dualDerivative, dualError, 1e-8);
	double time = 0;
	print(simulationDescriptor, " : entering simulation");
	const static double totalSimTime = 24.0 * 365 * 5 * 3600;
	clock_t time0 = clock();

	int previousPercent = 0;
	while (time < totalSimTime)
	{
		double percentComplete = time / totalSimTime;
		if ((int)(percentComplete * 100) > previousPercent)
		{
			previousPercent = (int)(percentComplete * 100);
			double ellapsedTime = (double)(clock() - time0) / CLOCKS_PER_SEC;
			double avgRate = percentComplete / ellapsedTime;
			double eta = (1 - percentComplete) / avgRate;
			print(simulationDescriptor, " : Percent Complete: ", (int)(percentComplete * 100));
			print(simulationDescriptor, " : ETA: ", formatTime(eta));
			print();
		}
		time += simulation.nextState();
	}
	print(simulationDescriptor, " : Simulation Complete");
	double ellapsedTime = (double)(clock() - time0) / CLOCKS_PER_SEC;
	print(simulationDescriptor, " : Ellapsed Time: ", formatTime(ellapsedTime));


	print(simulationDescriptor, " : formatting data");
	vector<double> timeList;
	vector<LISA_Vector> stateData;
	vector<LISA_Vector> perturbationData;
	for (totalState item : simulation.getData())
	{
		timeList.push_back(project<0>(item));
		stateData.push_back(project<1>(item));
		perturbationData.push_back(project<2>(item) - project<1>(item));
	}

	vector<vector<double>> finalData;
	for (int i = 0; i < stateData.size(); i++)
	{
		vector<double> temp1({timeList[i] / 365 / 24 / 3600 + 2029});
		LISA_Vector state(stateData[i]);
		LISA_Vector perturbation(perturbationData[i]);
		auto temp2 = flattendLISA(state);
		temp1.insert(temp1.end(), temp2.begin(), temp2.end());
		auto temp3 = flattendLISA(perturbation);
		temp1.insert(temp1.end(), temp3.begin(), temp3.end());
		finalData.push_back(temp1);
	}
	writeCSVfile(finalData, headerStringList,
		"Secular_Output_Directory/" + simulationDescriptor + "__Output.csv");
}

struct simulationVariables
{
	inputFileHandler data;
	vector<int> asteroidIndices;
	double leadingAngle;
	string prefixName;
	double maxTimeOffset;
	double deltaTimeOffset;
	double startingEpoch;
	int deltaAngle;
};

// void runLaggingBatch(inputFileHandler& data, vector<int>& asteroidIndices)
void runBatch(simulationVariables& vars)
{
	const static double earthOrbitalParameterEpoch = 28324.75;

	vector<vector<double>> anglesFromCeres;
	int timeIndex = 0;

	print("Starting Simulations for Lagging Batch");

	for (double timeOffset = 0;
		timeOffset < vars.maxTimeOffset;
		timeOffset += vars.deltaTimeOffset, ++timeIndex)
	{
		//set up initial values for everything but LISA itself. This
		//includes the asteroids and earth. This has to be done
		//inside this loop because these depend on the initial time
		//of the simulation.
		double startingTime = vars.startingEpoch * 24 * 3600 + timeOffset;
		asteroids = vars.data.generateOrbits(vars.asteroidIndices, startingTime);

		asteroidPositions = vector<asteroidData>();
		for (auto&& asteroid : asteroids)
		{
			asteroidPositions.push_back(asteroidData(asteroid.pos(), asteroid.getMass()));
		}

		Orbit earth(radians(200.7),//mean anomaly
			0.01671,//eccentricity
			AU,//semi major axis
			radians(0),//inclination
			radians(-11.261),//longitude of ascending node
			radians(114.2078),//argument of  the perihelion
			true,//provided mean anomaly, not true anomaly
			earthOrbitalParameterEpoch * 24 * 3600,//parameter epoch
			startingTime,//initial simulation epoch
			5.972e24,//Earth mass
			"Earth");//name

		double longitudeOfLISA = atan2(earth.y(), earth.x()) + PI * vars.leadingAngle / 180;
		double longitudeOfCeres = atan2(asteroids[0].y(), asteroids[0].x());
		anglesFromCeres.push_back({
			longitudeOfLISA - longitudeOfCeres
		});

		for (double constellationAngle = 0; constellationAngle < 120; constellationAngle += vars.deltaAngle)
		{
			LISA lisa(constellationAngle * PI / 180, longitudeOfLISA);

			runSimulation(lisa, vars.prefixName + "__Angle_" + to_string((int)constellationAngle)
				+ "__Time_Index_" + to_string(timeIndex));
		}
	}

	writeCSVfile(anglesFromCeres,
		{"angles_from_ceres"},
		"Secular_Output_Directory/" + vars.prefixName + "__Sampled_Angles_from_Ceres.csv");

	vector<vector<double>> constellationAngles;
	for (double constellationAngle = 0; constellationAngle < 120; constellationAngle += vars.deltaAngle)
	{
		constellationAngles.push_back({constellationAngle});
	}
	writeCSVfile(constellationAngles,
		{"initial_constellation_angles"},
		"Secular_Output_Directory/" + vars.prefixName + "__Sampled_Constellation_Angles.csv");
}

int main(int argc, char** argv){
	print("Initializing");
	inputFileHandler data(1000000);//start with 1000000 gaussian numbers
	vector<int> asteroidIndices(1);//we will use 1 asteroid
	for (int i = 0; i < asteroidIndices.size(); ++i)
	{
		asteroidIndices[i] = i;//fill array of indices of the asteroids we want
	}


	double fullRotationTime = 466.6 * 24 * 3600;
	double maxTimeOffset = 3 * fullRotationTime;
	int deltaAngle = 5;
	// double maxTimeOffset = fullRotationTime;
	// int deltaAngle = 20;
	double deltaTimeOffset = fullRotationTime * deltaAngle / 360;
	double startingEpoch = yearToEpoch(2029);

	simulationVariables vars{
		data,
		asteroidIndices,
		-20,
		"Lagging",
		maxTimeOffset,
		deltaTimeOffset,
		startingEpoch,
		deltaAngle
	};

	runBatch(vars);
	vars.leadingAngle = 20;
	vars.prefixName = "Leading";
	runBatch(vars);

	return 0;
}


