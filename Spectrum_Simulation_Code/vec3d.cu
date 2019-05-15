#include "vec3d.cuh"

#include <string>

#include "globalFunctions.h"

#define TOTAL_ACCESS __host__ __device__

using namespace std;
using namespace cudaUtil;

TOTAL_ACCESS vec3d::vec3d()
{
	x = 0;
	y = 0;
	z = 0;
}
TOTAL_ACCESS vec3d::vec3d(float inX, float inY)
{
	x = inX;
	y = inY;
	z = 0;
}
TOTAL_ACCESS vec3d::vec3d(float inX, float inY, float inZ)
{
	x = inX;
	y = inY;
	z = inZ;
}
TOTAL_ACCESS vec3d & vec3d::operator+=(const vec3d rhs)
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
	return *this;
}
TOTAL_ACCESS vec3d & vec3d::operator-=(const vec3d rhs)
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
	return *this;
}
TOTAL_ACCESS vec3d & vec3d::operator*=(const float rhs)
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
	return *this;
}
TOTAL_ACCESS vec3d & vec3d::operator/=(const float & rhs)
{
	x /= rhs;
	y /= rhs;
	z /= rhs;
	return *this;
}
string vec3d::toString()
{
	return "< " + str(x) + "," + str(y)
		+ "," + str(z) + " >";
}
TOTAL_ACCESS float vec3d::dot(vec3d otherVector)
{
	return x * otherVector.x + y * otherVector.y + z * otherVector.z;
}
TOTAL_ACCESS float vec3d::magnitude()
{
	return sqrt(x*x + y*y + z*z);
}
TOTAL_ACCESS float vec3d::squared()
{
	return x*x + y*y + z*z;
}
TOTAL_ACCESS vec3d vec3d::normalized()
{
	return vec3d(*this) / magnitude();
}
vector<float> vec3d::toFloatVector()
{
	vector<float> output(3);
	output[0] = x;
	output[1] = y;
	output[2] = z;
	return output;
}

TOTAL_ACCESS vec3d cudaUtil::operator+(vec3d lhs, const vec3d & rhs)
{
	lhs += rhs;
	return lhs;
}
TOTAL_ACCESS vec3d cudaUtil::operator-(vec3d lhs, const vec3d & rhs)
{
	lhs -= rhs;
	return lhs;
}
TOTAL_ACCESS vec3d cudaUtil::operator*(vec3d lhs, const float & rhs)
{
	lhs *= rhs;
	return lhs;
}
TOTAL_ACCESS vec3d cudaUtil::operator*(const float & lhs, vec3d rhs)
{
	rhs *= lhs;
	return rhs;
}
TOTAL_ACCESS vec3d cudaUtil::operator/(vec3d lhs, const float & rhs)
{
	lhs /= rhs;
	return lhs;
}
