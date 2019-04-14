#pragma once

#include <string>
#include "VectorSpace.h"

namespace lisa{
	class Orbit
	{
		double trueAnomaly, eccentricAnomaly, meanAnomaly, meanAnomaly0;
		double semiMajorAxis, semiLatusRectum, eccentricity;
		double inclination, longitudeOfAscendingNode, argumentOfPeriapsis;
		double meanAngularMotion, time, mass;
		double cLongAsc, sLongAsc, cInc, sInc;
		double r;
		std::string name;
		inline void trueToEccAnom();
		inline void eccToTrueAnom();
		inline void meanToEccAnom();
		inline void makeR();
		inline void eccToMeanAnom();
	public:
		Orbit();
		Orbit(double anomaly, double ecc, double semiMajor,
			double incline, double longAscend,
			double argPeri, bool isMeanAnom);
		Orbit(double anomaly, double ecc, double semiMajor,
			double incline, double longAscend,
			double argPeri, bool isMeanAnom, double timeOfParameters,
			double definedTimeZero, double inMass, std::string inName);
		void setTime(double newTime);
		double x();
		double y();
		double z();
		vsu::ProductSpace<double, double, double> pos();
		double vx();
		double vy();
		double vz();
		vsu::ProductSpace<double, double, double> vel();
		double getMass();
		std::string getName();
	};
}

#undef TOTAL_ACCESS
