import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import my_package.util as tools
import math

from Simulation_Core import *
from Data_Extraction_Methods import *

simulationDataDirectory = "../../Build_Directories/LISA/Secular_Output_Directory/"

prefixes = ["Lagging", "Leading"]

def graphMesh(xList, yList, inputMesh, graphTitle):
	levels = list(range(5, 35, 5))
	fmt = {}
	for level in levels:
		fmt[level] = "%.1f"%(level)
	CS = plt.contour(xList, yList, inputMesh, levels = levels)
	plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)

	plt.title(graphTitle)
	plt.xlabel("Initial Angle to Ceres")
	plt.ylabel("Initial Constellation Angle")
	plt.show()


def findMaxPerturbation(stacks):
	armDifferences = getD_ArmLengths(stacks)
	maxDifs = []
	for singleArmDifferences in armDifferences:
		maxDifs.append(0)
		for singleDifference in singleArmDifferences:
			maxDifs[-1] = max(maxDifs[-1], abs(singleDifference))
	return maxDifs

def stackMeshes(meshList):
	totalMesh = []
	for mesh in meshList:
		totalMesh.extend(mesh)
	return totalMesh

for prefix in prefixes:
	ceresAngleLocation = simulationDataDirectory\
		+ prefix + "__Sampled_Angles_from_Ceres.csv"
	constellationAngleLocation = simulationDataDirectory\
		+ prefix + "__Sampled_Constellation_Angles.csv"

	constellationAngles = tools.transpose(
		tools.floatParseCSV(constellationAngleLocation))[0]
	ceresAngles = tools.transpose(
		tools.floatParseCSV(ceresAngleLocation))[0]
	ceresAngleIndices = list(range(len(ceresAngles)))

	perturbationMeshes = [[],[],[]]
	for constellationAngle in constellationAngles:
		perturbationLists = [[],[],[]]
		for ceresAngleIndex in ceresAngleIndices:
			fileLocation = simulationDataDirectory\
				+ prefix + "__Angle_"\
				+ str(int(constellationAngle)) + "__Time_Index_"\
				+ str(ceresAngleIndex) + "__Output.csv"
			(stacks, times) = buildStacksAndTimes(fileLocation)

			perturbations = findMaxPerturbation(stacks)
			perturbationLists[0].append(perturbations[0])
			perturbationLists[1].append(perturbations[1])
			perturbationLists[2].append(perturbations[2])
		perturbationMeshes[0].append(perturbationLists[0])
		perturbationMeshes[1].append(perturbationLists[1])
		perturbationMeshes[2].append(perturbationLists[2])

	for i in range(len(ceresAngles)):
		ceresAngles[i] *= 180 / np.pi

	for i in range(1, len(ceresAngles)):
		while abs(ceresAngles[i] - ceresAngles[i-1]) > 180:
			ceresAngles[i] += 360

	bigMesh = stackMeshes(perturbationMeshes)
	graphMesh(ceresAngles, list(range(0, 360, 20)), bigMesh,
		"Maximum Perturbation Experienced on Arbitrary Arm "
			+ prefix + " Earth.")
	for i in range(3):
		graphMesh(ceresAngles, constellationAngles,
			perturbationMeshes[i],
			"Maximum Perturbation Experienced on Arm " + str(i + 1)
				+ " " + prefix + " Earth.")
	
