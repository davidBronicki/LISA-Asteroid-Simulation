#pragma once

#include <string>
#include <vector>

#include "cuda_runtime.h"

#define TOTAL_ACCESS __host__ __device__

namespace cudaUtil{
	class vec3d
	{
		float x;
		float y;
		float z;
	public:
		TOTAL_ACCESS vec3d();
		TOTAL_ACCESS vec3d(float inX, float inY);
		TOTAL_ACCESS vec3d(float inX, float inY, float inZ);
		TOTAL_ACCESS vec3d & operator+=(const vec3d rhs);
		TOTAL_ACCESS friend vec3d operator+(vec3d lhs, const vec3d & rhs);
		TOTAL_ACCESS vec3d & operator-=(const vec3d rhs);
		TOTAL_ACCESS friend vec3d operator-(vec3d lhs, const vec3d & rhs);
		TOTAL_ACCESS vec3d & operator*=(const float rhs);
		TOTAL_ACCESS friend vec3d operator*(vec3d lhs, const float & rhs);
		TOTAL_ACCESS friend vec3d operator*(const float & lhs, vec3d rhs);
		TOTAL_ACCESS vec3d & operator/=(const float & rhs);
		TOTAL_ACCESS friend vec3d operator/(vec3d lhs, const float & rhs);
		std::string toString();
		TOTAL_ACCESS float dot(vec3d otherVector);
		TOTAL_ACCESS float magnitude();
		TOTAL_ACCESS float squared();
		TOTAL_ACCESS vec3d normalized();
		std::vector<float> toFloatVector();
	};
	TOTAL_ACCESS vec3d operator+(vec3d lhs, const vec3d & rhs);
	TOTAL_ACCESS vec3d operator-(vec3d lhs, const vec3d & rhs);
	TOTAL_ACCESS vec3d operator*(vec3d lhs, const float & rhs);
	TOTAL_ACCESS vec3d operator*(const float & lhs, vec3d rhs);
	TOTAL_ACCESS vec3d operator/(vec3d lhs, const float & rhs);
}

#undef TOTAL_ACCESS
