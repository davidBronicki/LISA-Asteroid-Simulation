#include "LISA.h"

#include <vector>
#include <iostream>

#include "Orbit.h"
#include "VectorSpace.h"

using namespace std;
using namespace lisa;
using namespace vsu;

#define print(input) cout << input << endl

#include "defineOrbitalValues.h"

typedef ProductSpace<double, double, double> vec3;

LISA::LISA(double initialOrientiationAngle,
		double initialLongitude){
	static double alpha = ARM_LENGTH / (2 * AU);
	static double x = 2/sqrt(3) * alpha;
	static double ecc = sqrt(1 + x + x * x) - 1;
	static double incline = atan(alpha / (1 + alpha / sqrt(3)));
	static double orbitalOffsetAngle = TWOPI / 3;
	for (int i = 0; i < 3; ++i)
	{
		satellites.push_back(Orbit(
			HALFPI + initialOrientiationAngle + i * orbitalOffsetAngle + initialLongitude,//mean anomaly
			ecc,//eccentricity
			AU,//semi major axis
			incline,//inclination
			-initialOrientiationAngle - i * orbitalOffsetAngle,//longitude of ascending node
			-HALFPI,//argument of perihelion
			true, 0, 0, 1, "LISA Craft " + to_string(i + 1)));
	}
}

LISA::LISA(double initialOrientiationAngle)
:
	LISA(initialOrientiationAngle, -20.0/180 * PI)
{}

LISA::LISA()
:
	LISA(0)
{}

void LISA::setTime(double time){
	satellites[0].setTime(time);
	satellites[1].setTime(time);
	satellites[2].setTime(time);
}

vector<vec3> LISA::getPositions(){
	return vector<vec3>({satellites[0].pos(), satellites[1].pos(), satellites[2].pos()});
}

vector<Orbit> LISA::getSats(){
	return satellites;
}

vector<vec3> LISA::getUnitVectors(vector<vec3> lisaPositions){
	vec3 disp12(lisaPositions[1] - lisaPositions[0]);
	vec3 disp13(lisaPositions[2] - lisaPositions[0]);
	vec3 disp23(lisaPositions[2] - lisaPositions[1]);
	return vector<vec3>({disp12 / disp12.defaultMagnitude(),
		disp13 / disp13.defaultMagnitude(),
		disp23 / disp23.defaultMagnitude()});
}

vector<vec3> LISA::getUnitVectors(){
	return getUnitVectors(getPositions());
}
