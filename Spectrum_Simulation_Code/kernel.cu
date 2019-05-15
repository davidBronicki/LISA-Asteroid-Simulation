#include "globalSwitchStatements.h"
#include "kernel.cuh"

#include "cuda_runtime.h"

// #include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec3d.cuh"
#include "Orbit.cuh"
#include "LISA.h"
#include "dataHandler.h"
#include "hdmem.cuh"

using namespace std;
using namespace cudaUtil;
using namespace lisa;

#ifdef DEBUG_MODE
#define print(input) cout << "kernel: " << input << endl
#else
#ifdef NO_PRINT
#define print(input)
#else
#define print(input) cout << input << endl
#endif
#endif

#define syncThreads __syncthreads()
#define syncDevice cudaDeviceSynchronize()

#define THREAD_COUNT 128

#define G 6.67408e-11f

__global__ void coreKernel(Orbit* asteroids,
	vec3d* lisaLocs, vec3d* lisaUnitVectors,
	float time, vec3d* forceRelOutput, vec3d* forceAbsOutput,
	int n, int asteroidCount, int blocksPerSampleTime)
{
	//by sharing memory, a summation can be performed
	//within this method call to save on memory
	__shared__ vec3d tempVals[THREAD_COUNT];
	__shared__ vec3d tempAbsVals[3 * THREAD_COUNT];
	syncThreads;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < (asteroidCount))
	{
		float mu = asteroids[index].getMass() * G;//M*G
		asteroids[index].setTime(time);//update asteroid to current time
		vec3d pos = asteroids[index].pos();//get updated position
		//get vector-wise displacement from each satellite
		vec3d disp1 = lisaLocs[0] - pos;
		vec3d disp2 = lisaLocs[1] - pos;
		vec3d disp3 = lisaLocs[2] - pos;
		//get scalor distance to each satellite
		float dist1 = disp1.magnitude();
		float dist2 = disp2.magnitude();
		float dist3 = disp3.magnitude();
		//get gravitational acceleration on each satellite
		tempAbsVals[3 * threadIdx.x] = -mu * disp1 / (dist1 * dist1 * dist1);
		tempAbsVals[3 * threadIdx.x + 1] = -mu * disp2 / (dist2 * dist2 * dist2);
		tempAbsVals[3 * threadIdx.x + 2] = -mu * disp3 / (dist3 * dist3 * dist3);
		//get projection of accelerations onto each arm
		tempVals[threadIdx.x] =
			vec3d((tempAbsVals[3 * threadIdx.x + 1] - tempAbsVals[3 * threadIdx.x]).dot(lisaUnitVectors[0]),
			(tempAbsVals[3 * threadIdx.x + 2] - tempAbsVals[3 * threadIdx.x]).dot(lisaUnitVectors[1]),
			(tempAbsVals[3 * threadIdx.x + 2] - tempAbsVals[3 * threadIdx.x + 1]).dot(lisaUnitVectors[2]));
	}
	syncThreads;
	int i = blockDim.x;
	//sum over all threads in this block
	while (i != 1)
	{
		int j = i & 1;
		i >>= 1;
		if (threadIdx.x < i)
		{
			tempVals[threadIdx.x] += tempVals[threadIdx.x + i + j];
			tempAbsVals[3 * threadIdx.x] += tempAbsVals[3 * (threadIdx.x + i + j)];
			tempAbsVals[3 * threadIdx.x + 1] += tempAbsVals[3 * (threadIdx.x + i + j) + 1];
			tempAbsVals[3 * threadIdx.x + 2] += tempAbsVals[3 * (threadIdx.x + i + j) + 2];
		}
		i += j;
		syncThreads;
	}
	//set output value
	if (threadIdx.x == 0)
	{
		forceRelOutput[blockIdx.x + (n) * (blocksPerSampleTime)] = tempVals[0];
		forceAbsOutput[3 * (blockIdx.x + (n) * (blocksPerSampleTime))] = tempAbsVals[0];
		forceAbsOutput[3 * (blockIdx.x + (n) * (blocksPerSampleTime)) + 1] = tempAbsVals[1];
		forceAbsOutput[3 * (blockIdx.x + (n) * (blocksPerSampleTime)) + 2] = tempAbsVals[2];
	}
}

__global__ void compactingKernel(
	vec3d* forceRelOutput,
	vec3d* compactForceRelOutput,
	vec3d* forceAbsOutput,
	vec3d* compactForceAbsOutput,
	int sampleTimes, int blocksPerSampleTime)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < sampleTimes)
	{
		compactForceRelOutput[index] = vec3d();
		compactForceAbsOutput[3 * index] = vec3d();
		compactForceAbsOutput[3 * index + 1] = vec3d();
		compactForceAbsOutput[3 * index + 2] = vec3d();
		for (int i = 0; i < (blocksPerSampleTime); i++)
		{
			compactForceRelOutput[index] +=
				forceRelOutput[index * (blocksPerSampleTime) + i];
			compactForceAbsOutput[3 * index] +=
				forceAbsOutput[3 * (index * (blocksPerSampleTime) + i)];
			compactForceAbsOutput[3 * index + 1] +=
				forceAbsOutput[3 * (index * (blocksPerSampleTime) + i) + 1];
			compactForceAbsOutput[3 * index + 2] +=
				forceAbsOutput[3 * (index * (blocksPerSampleTime) + i) + 2];
		}
	}
}

bool executeCudaCode(dataHandler& data)
{
	vector<float> sampleTimes = data.generateSampleTimes();
	vector<Orbit> asteroids = data.generateAsteroids();
	LISA lisa = data.generateLISA();
	//final output vectors
	vector<vec3d> outputRel(sampleTimes.size());
	vector<vector<vec3d>> outputAbs(sampleTimes.size());
	//vector<vector<vec3d>> unitVectors(sampleTimeCount);

	//number of blocks required to simulate a single intance of time
	int blocksPerSampleTime = (asteroids.size() + THREAD_COUNT - 1) / THREAD_COUNT;

	//three vectors for position and three vectors for unit vectors
	//for LISA at each point of time. First dimension specifies time,
	//second dimension specifies which of three.
	vector<vector<vec3d>> lisaLocationsVector(sampleTimes.size());
	vector<vector<vec3d>> lisaUnitVectorsVector(sampleTimes.size());

	for (int i = 0; i < sampleTimes.size(); i++){
		lisa.setTime(sampleTimes[i]);
		lisaLocationsVector[i] = lisa.getPositions();
		lisaUnitVectorsVector[i] = lisa.getUnitVectors(lisaLocationsVector[i]);
	}

	// //device memory instantiated with asteroids vector
	hdmem<Orbit> orbits = hdmem<Orbit>(asteroids).loadMemory().updateDevice();

	//single time slot of lisa data
	hdmem<vec3d> lisaLocations = hdmem<vec3d>(3).loadMemory();
	hdmem<vec3d> lisaUnitVectors = hdmem<vec3d>(3).loadMemory();
	//device memory of output stuff
	hdmem<vec3d> partialRelResults =
		hdmem<vec3d>(sampleTimes.size() * blocksPerSampleTime).loadMemory();
	hdmem<vec3d> partialAbsResults = 
		hdmem<vec3d>(3 * sampleTimes.size() * blocksPerSampleTime).loadMemory();
	hdmem<vec3d> relResults = hdmem<vec3d>(sampleTimes.size()).loadMemory();
	hdmem<vec3d> absResults = hdmem<vec3d>(3 * sampleTimes.size()).loadMemory();

	print((float)(orbits.deviceMemoryAllocated()
		+ lisaLocations.deviceMemoryAllocated()
		+ lisaUnitVectors.deviceMemoryAllocated()
		+ partialRelResults.deviceMemoryAllocated()
		+ partialAbsResults.deviceMemoryAllocated()
		+ relResults.deviceMemoryAllocated()
		+ absResults.deviceMemoryAllocated())/1000000 << " MB Allocated on Device");

	//Core simulation loop
	dim3 grid(blocksPerSampleTime);
	dim3 block(THREAD_COUNT);
	for (int i = 0; i < sampleTimes.size(); i++)
	{
		//update time specific values and check validity of device memory copy
		bool working = true;
		working &= lisaLocations.setAll(lisaLocationsVector[i]).isGood();
		working &= lisaUnitVectors.setAll(lisaUnitVectorsVector[i]).isGood();
		if (!working)
		{
			print("failure on step " << i << ". (0 indexed)");
			print("lisa location status");
			lisaLocations.updateHost();
			print("lisu unit vectors status");
			lisaUnitVectors.updateHost();
			print("orbit status");
			orbits.updateHost();
			print("partial relative results status");
			partialRelResults.updateHost();
			print("partial absolute results status");
			partialAbsResults.updateHost();
			print("relative results status");
			relResults.updateHost();
			print("absolute results status");
			absResults.updateHost();
			// killProgram();
			goto killSpace;
		}
		//simulate single point of time
		coreKernel << <grid, block >> >
			(orbits, lisaLocations, lisaUnitVectors, sampleTimes[i],
				partialRelResults, partialAbsResults,
				i, orbits.size(), blocksPerSampleTime);
		syncDevice;
	}
	//perform final sumation for each point of time.
	grid = dim3((sampleTimes.size() + THREAD_COUNT - 1) / THREAD_COUNT);
	block = dim3(THREAD_COUNT);
	compactingKernel << <grid, block >> >
		(partialRelResults, relResults, partialAbsResults,
			absResults, sampleTimes.size(), blocksPerSampleTime);
	syncDevice;
	//fill in output vectors
	relResults.updateHost();
	absResults.updateHost();
	for (int i = 0; i < sampleTimes.size(); i++)
	{
		outputRel[i] = relResults.get(i);
		for (int j = 0; j < 3; j++)
		{
			outputAbs[i].push_back(absResults.get(3 * i + j));
		}
	}
	data.handOffResults(outputRel, outputAbs, lisaUnitVectorsVector, true);

	return true;
killSpace:
	return false;
}