#pragma once

#include "cuda_runtime.h"
#include <string>
#include "vec3d.cuh"

#define TOTAL_ACCESS __host__ __device__

namespace lisa{
	class Orbit
	{
		float trueAnomaly;
		float eccentricAnomaly;
		float meanAnomaly;
		float semiMajorAxis;
		float semiLatusRectum;
		float eccentricity;
		float inclination;
		float longitudeOfAscendingNode;
		float argumentOfPeriapsis;
		float meanAngularMotion;
		float time;
		float time0;
		float meanAnomaly0;
		float mass;
		std::string name;
		float cLongAsc;
		float sLongAsc;
		float cInc;
		float sInc;
		float r;
		TOTAL_ACCESS void TrueToEccAnom();
		TOTAL_ACCESS void EccToTrueAnom();
		TOTAL_ACCESS void MeanToEccAnom();
	public:
		Orbit();
		Orbit(float anomaly, float ecc, float semiMajor,
			float incline, float longAscend,
			float argPeri, bool isTrueAnom);
		Orbit(float anomaly, float ecc, float semiMajor,
			float incline, float longAscend,
			float argPeri, bool isTrueAnom, float timeOfParameters,
			float definedTimeZero, float inMass, std::string inName);
		TOTAL_ACCESS void setTime(float newTime);
		TOTAL_ACCESS void setTime(float newTime, bool useAbsoluteTime);
		TOTAL_ACCESS float x();
		TOTAL_ACCESS float y();
		TOTAL_ACCESS float z();
		TOTAL_ACCESS float getMass();
		TOTAL_ACCESS cudaUtil::vec3d pos();
		std::string getName();
	};
}

#undef TOTAL_ACCESS
