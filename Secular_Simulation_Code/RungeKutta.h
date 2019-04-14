#pragma once

#include <vector>

template<typename T>
class RungeKuttaSimulation{
	T (*stateDerivative)(const T&);
	double (*numericalError)(const T&, const T&, double);
	std::vector<T> stateList;
	double errorParameter;//specifies how small of an error is acceptable,
	//allows quick manual change to error without changing numbericalError.
	//Value stored is inverse of parameter passed: the value passed is small
	//for high precision, value stored is high for high precision.
	double dt;
	inline T runStep(const T& inState, double timeStep)
	{
		auto s1 = stateDerivative(stateList.back());
		auto s2 = stateDerivative(stateList.back() + timeStep / 2 * s1);
		auto s3 = stateDerivative(stateList.back() + timeStep / 2 * s2);
		auto s4 = stateDerivative(stateList.back() + timeStep * s3);
		return inState + timeStep / 6 * ((s1 + 2*s2) + (2*s3 + s4));
	}
public:
	RungeKuttaSimulation(T initialState, double initialTimeStep,
		T (*inStateDerivative)(const T&),
		double (*inNumericalError)(const T&, const T&, double),
		double inErrorParameter)
	:
		dt(initialTimeStep),
		stateDerivative(inStateDerivative),
		numericalError(inNumericalError),
		errorParameter(1.0/inErrorParameter)
	{
		stateList.push_back(initialState);
	}
	double nextState(){
		auto testState = runStep(stateList.back(), dt * 2);
		stateList.push_back(runStep(stateList.back(), dt));
		stateList.push_back(runStep(stateList.back(), dt));
		testState -= stateList.back();
		double error = errorParameter * numericalError(stateList.back(), testState, dt);
		double timeEllapsed = 2 * dt;
		if (error > 1) dt *= pow(error, -.25);
		else dt *= pow(error, -.2);
		return timeEllapsed;
	}
	std::vector<T> getData(){
		return stateList;
	}
	T getState(){
		return stateList.back();
	}
};
