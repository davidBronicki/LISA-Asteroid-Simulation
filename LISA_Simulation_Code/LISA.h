#pragma once

#include <vector>
#include <math.h>

#include "Orbit.cuh"
#include "vec3d.cuh"

namespace lisa{
	class LISA{
		std::vector<Orbit> satellites;
	public:
		LISA();
		LISA(float initialOrientiationAngle);
		LISA(float initialOrientiationAngle,
			float initialLongitude);
		LISA(float initialOrientiationAngle,
			float initialLongitude,
			float initialTime);
		void setTime(float time);
		std::vector<cudaUtil::vec3d> getPositions();
		std::vector<cudaUtil::vec3d> getUnitVectors();
		std::vector<cudaUtil::vec3d> getUnitVectors(std::vector<cudaUtil::vec3d> lisPositions);
	};
}
