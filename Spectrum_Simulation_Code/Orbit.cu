#include "Orbit.cuh"

#include <math.h>
#include <iostream>

#include "vec3d.cuh"

using namespace std;
using namespace lisa;
using namespace cudaUtil;

#define makeR r = semiLatusRectum / (1 + eccentricity * cos(trueAnomaly))
#define totalAngle (trueAnomaly + argumentOfPeriapsis)
#define eccToMeanAnom meanAnomaly = eccentricAnomaly - eccentricity * sin(eccentricAnomaly)
#define difOfMean(psi) (meanAnomaly - psi + eccentricity * sin(psi))

#define print(input) cout << input << endl

#define TOLERANCE 1e-10f
#define PI 3.1415926535898f
#define TWOPI 6.2831853071796f
#define au 149597870700.0f
#define G 6.67408e-11f
#define earthMass 5.972e24f
#define solarMass 1.989e30f
#define mu 1.32747e20f
#define TOTAL_ACCESS __host__ __device__

TOTAL_ACCESS float floatMod(float value, float mod)
{
	float extra = value >= 0 ?
		mod * (int)(value / mod) :
		mod * (int)(value / mod - 1);
	return value - extra;
}

Orbit::Orbit()
{
	trueAnomaly = 0;
	eccentricAnomaly = 0;
	meanAnomaly = 0;
	meanAnomaly0 = 0;
	semiMajorAxis = au;
	semiLatusRectum = au;
	eccentricity = 0;
	inclination = 0;
	longitudeOfAscendingNode = 0;
	argumentOfPeriapsis = 0;
	meanAngularMotion = 1.991e-7f;
	cLongAsc = 1;
	sLongAsc = 0;
	cInc = 1;
	sInc = 1e-7f;
	r = au;

	time = 0;
	time0 = 0;
	mass = 1;
	name = "default name";
}

Orbit::Orbit(float anomaly, float ecc, float semiMajor,
	float incline, float longAscend,
	float argPeri, bool isTrueAnom)
{
	eccentricity = ecc;
	if (isTrueAnom)
	{
		trueAnomaly = floatMod(anomaly, TWOPI);
		TrueToEccAnom();
		eccToMeanAnom;
	}
	else
	{
		meanAnomaly = floatMod(anomaly, TWOPI);
		MeanToEccAnom();
		EccToTrueAnom();
	}
	meanAnomaly0 = meanAnomaly;
	semiMajorAxis = semiMajor;
	semiLatusRectum = semiMajor * (1 - ecc * ecc);
	inclination = incline;
	longitudeOfAscendingNode = floatMod(longAscend, TWOPI);
	argumentOfPeriapsis = floatMod(argPeri, TWOPI);
	meanAngularMotion = sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis));
	cLongAsc = cos(longitudeOfAscendingNode);
	sLongAsc = sin(longitudeOfAscendingNode);
	cInc = cos(inclination);
	sInc = sin(inclination);
	if (sInc == 0)
	{
		sInc = 1e-7;
	}
	makeR;

	mass = 1;
	name = "default name";
}
Orbit::Orbit(float anomaly, float ecc, float semiMajor,
	float incline, float longAscend,
	float argPeri, bool isTrueAnom, float timeOfParameters,
	float definedTimeZero, float inMass, std::string inName)
{
	eccentricity = ecc;
	if (isTrueAnom)
	{
		trueAnomaly = floatMod(anomaly, TWOPI);
		TrueToEccAnom();
		eccToMeanAnom;
	}
	else
	{
		meanAnomaly = floatMod(anomaly, TWOPI);
		MeanToEccAnom();
		EccToTrueAnom();
	}
	semiMajorAxis = semiMajor;
	semiLatusRectum = semiMajor * (1 - ecc * ecc);
	inclination = incline;
	longitudeOfAscendingNode = floatMod(longAscend, TWOPI);
	argumentOfPeriapsis = floatMod(argPeri, TWOPI);
	meanAngularMotion = sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis));
	cLongAsc = cos(longitudeOfAscendingNode);
	sLongAsc = sin(longitudeOfAscendingNode);
	cInc = cos(inclination);
	sInc = sin(inclination);
	if (sInc == 0)
	{
		sInc = 1e-7;
	}
	makeR;

	meanAnomaly0 = meanAnomaly + (definedTimeZero - timeOfParameters) * meanAngularMotion;
	meanAnomaly = meanAnomaly0;
	MeanToEccAnom();
	EccToTrueAnom();
	time0 = definedTimeZero;
	time = definedTimeZero;
	mass = inMass;
	name = inName;
}

TOTAL_ACCESS void Orbit::TrueToEccAnom()
{
	float temp = sqrt((1 - eccentricity) / (1 + eccentricity)) * tan(trueAnomaly / 2);
	temp = 2 * atan(temp);
	eccentricAnomaly = temp >= 0 ? temp : temp + TWOPI;
	//eccentricAnomaly = 2 * atan(temp);
}
TOTAL_ACCESS void Orbit::EccToTrueAnom()
{
	float temp = sqrt((1 + eccentricity) / (1 - eccentricity)) * tan(eccentricAnomaly / 2);
	temp = 2 * atan(temp);
	trueAnomaly = temp > 0 ? temp : temp + TWOPI;
	//trueAnomaly = 2 * atan(temp);
}
TOTAL_ACCESS void Orbit::MeanToEccAnom()
{
	//int temp = (int)(meanAnomaly / TWOPI);
	//meanAnomaly -= temp * TWOPI;
	double lower = 0;
	double upper = TWOPI;
	double width = TWOPI;
	double midPoint = PI;
	while (width > TOLERANCE)
	{
		if (difOfMean(midPoint) * difOfMean(lower) > 0)
		{
			lower = midPoint;
		}
		else
		{
			upper = midPoint;
		}
		width = upper - lower;
		midPoint = (upper + lower) / 2;
	}
	eccentricAnomaly = midPoint;
}

TOTAL_ACCESS void Orbit::setTime(float newTime)
{
	time = newTime;
	//meanAnomaly = meanAnomaly0 + (time - time0) * meanAngularMotion;
	meanAnomaly = floatMod(meanAnomaly0 + (time * meanAngularMotion), TWOPI);
	MeanToEccAnom();
	EccToTrueAnom();
	makeR;
}

TOTAL_ACCESS void Orbit::setTime(float newTime, bool useAbsoluteTime){
	if (useAbsoluteTime){
		time = newTime;
		//meanAnomaly = meanAnomaly0 + (time - time0) * meanAngularMotion;
		meanAnomaly = floatMod(meanAnomaly0 + ((time - time0) * meanAngularMotion), TWOPI);
		MeanToEccAnom();
		EccToTrueAnom();
		makeR;
	}
	else setTime(newTime);
}

TOTAL_ACCESS float Orbit::x()
{
	return r * (cLongAsc * cos(totalAngle) - cInc * sLongAsc * sin(totalAngle));
}

TOTAL_ACCESS float Orbit::y()
{
	return r * (sLongAsc * cos(totalAngle) + cInc * cLongAsc * sin(totalAngle));
}

TOTAL_ACCESS float Orbit::z()
{
	return r * (sInc * sin(totalAngle));
}

TOTAL_ACCESS float Orbit::getMass()
{
	return mass;
}

TOTAL_ACCESS vec3d Orbit::pos()
{
	return vec3d(x(), y(), z());
}

string Orbit::getName()
{
	return name;
}

