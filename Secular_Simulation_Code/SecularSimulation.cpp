#include "globalSwitchStatements.h"
#include <iostream>
#include <math.h>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>
#include <thread>

#include "VectorSpace.h"
#include "Orbit.h"
#include "globalFunctions.h"
#include "dataHandler.h"
#include "RungeKutta.h"

#include "defineOrbitalValues.h"

using namespace std;
using namespace vsu;
using namespace lisa;

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

typedef ProductSpace<double, double, double> vec3;

typedef ProductSpace<vec3, vec3> state;
typedef ProductSpace<state, state, state> LISA_Vector;
typedef ProductSpace<double, LISA_Vector, LISA_Vector> totalState;

typedef ProductSpace<vec3, double> asteroidData;

vector<Orbit> asteroids;
vector<asteroidData> asteroidPositions;

vec3 asteroidForce(const vec3& displacement, double mass){
	double distance = displacement.defaultMagnitude();
	return displacement * mass * G / pow(distance, 3);
}

vec3 centralForce(const vec3& position)
{
	return -1 * position * MU / pow(position.defaultMagnitude(), 3);
}

state oribitalDerivativeWithAsteroids(const state& input, double time){
	auto position = Projection<0>::get(input);
	auto velocity = Projection<1>::get(input);
	vec3 force;
	for (asteroidData asteroid : asteroidPositions)
	{
		force += asteroidForce(Projection<0>::get(asteroid) - position, Projection<1>::get(asteroid));
	}
	return state(velocity, force + centralForce(position));
}

state oribitalDerivativeWithoutAsteroids(const state& input, double time){
	auto position = Projection<0>::get(input);
	auto velocity = Projection<1>::get(input);
	return state(velocity, centralForce(position));
}

LISA_Vector LISA_DerivativeWithAsteroids(const LISA_Vector& input, double time){
	return LISA_Vector(oribitalDerivativeWithAsteroids(Projection<0>::get(input), time),
		oribitalDerivativeWithAsteroids(Projection<1>::get(input), time),
		oribitalDerivativeWithAsteroids(Projection<2>::get(input), time));
}

LISA_Vector LISA_DerivativeWithoutAsteroids(const LISA_Vector& input, double time){
	return LISA_Vector(oribitalDerivativeWithoutAsteroids(Projection<0>::get(input), time),
		oribitalDerivativeWithoutAsteroids(Projection<1>::get(input), time),
		oribitalDerivativeWithoutAsteroids(Projection<2>::get(input), time));
}

#define THREAD_COUNT 2

totalState dualDerivative(const totalState& input){
	double time = Projection<0>::get(input);
	vector<thread> threads;
	for (int i = 0; i < THREAD_COUNT; i++){
		threads.push_back(thread([](int threadNumber, double time){
				for (int i = threadNumber; i < asteroids.size(); i += THREAD_COUNT){
					asteroids[i].setTime(time);
					asteroidPositions[i] = asteroidData(asteroids[i].pos(), asteroids[i].getMass());
				}
			}, i , time));
	}
	for (int i = 0; i < THREAD_COUNT; i++){
		threads[i].join();
	}
	return totalState(1,
		LISA_DerivativeWithoutAsteroids(Projection<1>::get(input), time),
		LISA_DerivativeWithAsteroids(Projection<2>::get(input), time));
}

double stateSqrError(const vec3& currentState, const vec3& deltaState, double dt)
{
	return deltaState.defaultSquareMagnitude() / currentState.defaultSquareMagnitude();
}

double orbitalSqrError(const state& currentState, const state& deltaState, double dt)
{
	double e1 = stateSqrError(Projection<0>::get(currentState), Projection<0>::get(deltaState), dt);
	double e2 = stateSqrError(Projection<1>::get(currentState), Projection<1>::get(deltaState), dt);
	return e1 + e2;
}

double LISA_SqrError(const LISA_Vector& currentState, const LISA_Vector& deltaState, double dt)
{
	double e1 = orbitalSqrError(Projection<0>::get(currentState), Projection<0>::get(deltaState), dt);
	double e2 = orbitalSqrError(Projection<1>::get(currentState), Projection<1>::get(deltaState), dt);
	double e3 = orbitalSqrError(Projection<2>::get(currentState), Projection<2>::get(deltaState), dt);
	return e1 + e2 + e3;
}

double dualError(const totalState& currentState, const totalState& deltaState, double dt)
{
	double e1 = LISA_SqrError(Projection<1>::get(currentState), Projection<1>::get(deltaState), dt);
	double e2 = LISA_SqrError(Projection<2>::get(currentState), Projection<2>::get(deltaState), dt);
	return sqrt(e1 + e2);
}


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

vector<double> flattendLISA(LISA_Vector input)
{
	vector<double> output;
	vector<state> sats({Projection<0>::get(input),
		Projection<1>::get(input),
		Projection<2>::get(input)});
	vector<vec3> vecs;
	for (auto sat : sats)
	{
		vecs.push_back(Projection<0>::get(sat));
		vecs.push_back(Projection<1>::get(sat));
	}
	for (auto vec : vecs)
	{
		output.push_back(Projection<0>::get(vec));
		output.push_back(Projection<1>::get(vec));
		output.push_back(Projection<2>::get(vec));
	}
	return output;
}

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

int main(int argc, char** argv){
	print("initializing");
	dataHandler data(argc, argv);
	asteroids = data.generateAsteroids();
	for (auto asteroid : asteroids)
	{
		// print("initial asteroid data: ", asteroidData(asteroid.pos(), asteroid.getMass()));
		asteroidPositions.push_back(asteroidData(asteroid.pos(), asteroid.getMass()));
	}
	print("Number of Asteroids: ");
	print(asteroidPositions.size());
	auto item = data.generateLISA();
	RungeKuttaSimulation<totalState> simulation(extractLISA(item), 1.0e4, dualDerivative, dualError, 1e-6);
	// RungeKuttaSimulation<totalState> simulation(extractLISA(item), 1.0e4, dualDerivative, dualError, 1e-10);
	double time = 0;
	print("entering simulation");
	// simulation.nextState();
	double totalSimTime = 24.0 * 365 * 5 * 3600;
	clock_t time0 = clock();
	while (time < totalSimTime)
	{
		double percentComplete = time / totalSimTime;
		double ellapsedTime = (double)(clock() - time0) / CLOCKS_PER_SEC;
		double avgRate = percentComplete / ellapsedTime;
		double eta = (1 - percentComplete) / avgRate;
		print("Percent Complete: ", percentComplete * 100);
		print("ETA: ", formatTime(eta));
		print();
		time += simulation.nextState();
	}
	print("Simulation Complete");
	double ellapsedTime = (double)(clock() - time0) / CLOCKS_PER_SEC;
	print("Ellapsed Time: ", formatTime(ellapsedTime));
	print("formatting data");
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
	writeCSVfile(finalData, "outputData.csv");
	return 0;
}

// 2029.000000

// 0006426875858,0149065164777,-666067850
// -29828,001409,00208

// 0125880781586,-080098420149,-666067850
// 016134,025128,00208

// -132307657445,-068966744627,-666067850
// 013694,-26537,00208
