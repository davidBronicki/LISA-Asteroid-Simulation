#include "LISA.h"

#include <vector>
#include <iostream>

#include "Orbit.cuh"
#include "vec3d.cuh"

using namespace std;
using namespace lisa;
using namespace cudaUtil;

#define print(input) cout << input << endl

#include "defineOrbitalValues.h"

LISA::LISA(float initialOrientiationAngle,
		float initialLongitude,
		float initialTime){
	satellites = vector<Orbit>(3);
	float alpha = ARM_LENGTH / (2 * AU);
	float ecc = sqrt(1 + (2 / sqrt(3)) * alpha + (4 / 3) * alpha * alpha) - 1;
	float incline = atan(alpha / (1 + alpha / (float)sqrt(3)));
	float satAngle = TWOPI / 3;
	float baseAngle = initialLongitude * TWOPI / 360;
	satellites[0] = Orbit(baseAngle + HALFPI,
		ecc, AU, incline, 0, -HALFPI,
		false, initialTime, initialTime, 1, "LISA Craft 1");
	satellites[1] = Orbit(baseAngle + HALFPI - satAngle,
		ecc, AU, incline, satAngle, -HALFPI,
		false, initialTime, initialTime, 1, "LISA Craft 2");
	satellites[2] = Orbit(baseAngle + HALFPI - 2 * satAngle,
		ecc, AU, incline, 2 * satAngle, -HALFPI,
		false, initialTime, initialTime, 1, "LISA Craft 3");
}

LISA::LISA(float initialOrientiationAngle,
		float initialLongitude):
		LISA(initialOrientiationAngle, initialLongitude, 0){}

LISA::LISA(float initialOrientiationAngle):
		LISA(initialOrientiationAngle, -20){}

LISA::LISA():LISA(0){}

void LISA::setTime(float time){
	satellites[0].setTime(time);
	satellites[1].setTime(time);
	satellites[2].setTime(time);
}

vector<vec3d> LISA::getPositions(){
	return vector<vec3d>({satellites[0].pos(), satellites[1].pos(), satellites[2].pos()});
}

vector<vec3d> LISA::getUnitVectors(vector<vec3d> lisaPositions){
	vec3d disp12 = lisaPositions[1] - lisaPositions[0];
	vec3d disp13 = lisaPositions[2] - lisaPositions[0];
	vec3d disp23 = lisaPositions[2] - lisaPositions[1];
	return vector<vec3d>({disp12.normalized(), disp13.normalized(), disp23.normalized()});
}

vector<vec3d> LISA::getUnitVectors(){
	return getUnitVectors(getPositions());
}
