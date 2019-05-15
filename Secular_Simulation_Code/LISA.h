#pragma once

#include <vector>
#include <math.h>

#include "Orbit.h"
#include "VectorSpace.h"

namespace lisa{
	class LISA{
		std::vector<Orbit> satellites;
	public:
		LISA();
		LISA(double initialOrientiationAngle);
		LISA(double initialOrientiationAngle,
			double initialLongitude);
		void setTime(double time);
		std::vector<vsu::ProductSpace<double, double, double>> getPositions();
		std::vector<Orbit> getSats();
		std::vector<vsu::ProductSpace<double, double, double>> getUnitVectors();
		std::vector<vsu::ProductSpace<double, double, double>> getUnitVectors(
			std::vector<vsu::ProductSpace<double, double, double>> lisaPositions);
	};
}
