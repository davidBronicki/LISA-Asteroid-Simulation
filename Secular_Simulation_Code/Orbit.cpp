#include "Orbit.h"

#include <math.h>
#include <iostream>

#include "VectorSpace.h"

using namespace std;
using namespace lisa;

#define totalAngle (trueAnomaly + argumentOfPeriapsis)
#define difOfMean(psi) (meanAnomaly - psi + eccentricity * sin(psi))

#define print(input) cout << input << endl

#define TOLERANCE 1e-10f

#include "defineOrbitalValues.h"

typedef vsu::ProductSpace<double, double, double> vec3;

inline double doubleMod(double value, double mod)
{
	double extra = value >= 0 ?
		mod * (int)(value / mod) :
		mod * (int)(value / mod - 1);
	return value - extra;
}

Orbit::Orbit()
:
	Orbit(0,0,AU,0,0,0,true)
{}

Orbit::Orbit(double anomaly, double ecc, double semiMajor,
	double incline, double longAscend,
	double argPeri, bool isTrueAnom)
:
	Orbit(anomaly, ecc, semiMajor, incline, longAscend,
		argPeri, isTrueAnom, 0, 0, 1, std::string("Default Name"))
{}

Orbit::Orbit(double anomaly, double ecc, double semiMajor,
	double incline, double longAscend,
	double argPeri, bool isMeanAnom, double timeOfParameters,
	double definedTimeZero, double inMass, std::string inName)
:
	eccentricity(ecc),
	semiMajorAxis(semiMajor),
	semiLatusRectum(semiMajor * (1-ecc*ecc)),
	inclination(incline),
	longitudeOfAscendingNode(doubleMod(longAscend, TWOPI)),
	argumentOfPeriapsis(doubleMod(argPeri, TWOPI)),
	meanAngularMotion(sqrt(MU/(semiMajor*semiMajor*semiMajor))),
	cLongAsc(cos(longAscend)),
	sLongAsc(sin(longAscend)),
	cInc(cos(incline)),
	sInc(sin(incline)),
	time(0),
	mass(inMass),
	name(inName)
{
	if (isMeanAnom)
	{
		meanAnomaly = doubleMod(anomaly, TWOPI);
	}
	else
	{
		trueAnomaly = doubleMod(anomaly, TWOPI);
		trueToEccAnom();
		eccToMeanAnom();
	}
	if (sInc == 0)
	{
		sInc = 1e-7;
	}

	meanAnomaly0 = doubleMod(meanAnomaly + (definedTimeZero - timeOfParameters)
		* meanAngularMotion, TWOPI);
	meanAnomaly = meanAnomaly0;
	meanToEccAnom();
	eccToTrueAnom();

	makeR();
}

inline void Orbit::trueToEccAnom()
{
	double temp = sqrt((1 - eccentricity) / (1 + eccentricity)) * tan(trueAnomaly / 2);
	temp = 2 * atan(temp);
	eccentricAnomaly = temp >= 0 ? temp : temp + TWOPI;
}
inline void Orbit::eccToTrueAnom()
{
	double temp = sqrt((1 + eccentricity) / (1 - eccentricity)) * tan(eccentricAnomaly / 2);
	temp = 2 * atan(temp);
	trueAnomaly = temp >= 0 ? temp : temp + TWOPI;
}
inline void Orbit::meanToEccAnom()
{
	eccentricAnomaly = meanAnomaly;
	double dif = TOLERANCE + 1;
	while (dif > TOLERANCE)
	{
		double prevEccAnom = eccentricAnomaly;
		eccentricAnomaly = meanAnomaly + eccentricity * sin(prevEccAnom);
		dif = abs(eccentricAnomaly - prevEccAnom);
	}
}

inline void Orbit::makeR()
{
	r = semiLatusRectum / (1 + eccentricity * cos(trueAnomaly));
}

inline void Orbit::eccToMeanAnom()
{
	meanAnomaly = eccentricAnomaly - eccentricity * sin(eccentricAnomaly);
}

void Orbit::setTime(double newTime)
{
	time = newTime;
	meanAnomaly = doubleMod(meanAnomaly0 + (time * meanAngularMotion), TWOPI);
	meanToEccAnom();
	eccToTrueAnom();
	makeR();
}

double Orbit::x()
{
	return r * (cLongAsc * cos(totalAngle) - cInc * sLongAsc * sin(totalAngle));
}

double Orbit::y()
{
	return r * (sLongAsc * cos(totalAngle) + cInc * cLongAsc * sin(totalAngle));
}

double Orbit::z()
{
	return r * (sInc * sin(totalAngle));
}

vec3 Orbit::pos()
{
	return vec3(x(), y(), z());
}

double Orbit::vx()
{
	return -sqrt(MU/semiLatusRectum) * (cLongAsc*(sin(argumentOfPeriapsis + trueAnomaly)
			+ eccentricity * sin(argumentOfPeriapsis))
		+ cInc * sLongAsc * (cos(argumentOfPeriapsis + trueAnomaly)
			+ eccentricity * cos(argumentOfPeriapsis)));
}

double Orbit::vy()
{
	return -sqrt(MU/semiLatusRectum) * (sLongAsc*(sin(argumentOfPeriapsis + trueAnomaly)
			+ eccentricity * sin(argumentOfPeriapsis))
		- cInc * cLongAsc * (cos(argumentOfPeriapsis + trueAnomaly)
			+ eccentricity * cos(argumentOfPeriapsis)));
}

double Orbit::vz()
{
	return sqrt(MU/semiLatusRectum) * sInc * (cos(argumentOfPeriapsis + trueAnomaly)
		+ eccentricity * cos(argumentOfPeriapsis));
}

vec3 Orbit::vel()
{
	return vec3(vx(), vy(), vz());
}

double Orbit::getMass()
{
	return mass;
}

string Orbit::getName()
{
	return name;
}

#include "undefineOrbitalValues.h"
