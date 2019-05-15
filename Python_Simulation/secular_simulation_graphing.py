import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import my_package.util as tools
import math

from Simulation_Core import *
from Data_Extraction_Methods import *

simulationDataDirectory = "../../Build_Directories/LISA/Secular_Output_Directory/"

prefixes = ["Lagging__", "Leading__"]

midWord = "__Time_Index_"

postfix = "__Output.csv"

constellationAngles = range(0, 120, 20)
ceresAngleIndices = range(0, 18)


def graphMesh(xList, yList, inputMesh, graphTitle):
	levels = list(range(5, 30, 5))
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

for prefix in prefixes:
	perturbationMeshes = [[],[],[]]
	for constellationAngle in constellationAngles:
		perturbationLists = [[],[],[]]
		for ceresAngleIndex in ceresAngleIndices:
			fileLocation = simulationDataDirectory\
				+ prefix + "Angle_"\
				+ str(constellationAngle) + midWord\
				+ str(ceresAngleIndex) + postfix
			(stacks, times) = buildStacksAndTimes(fileLocation)

			# print(fileLocation)
			# temp = getD_ArmLengths(stacks)
			# plt.plot(times, temp[0])
			# plt.plot(times, temp[1])
			# plt.plot(times, temp[2])
			# plt.grid()
			# plt.show()

			perturbations = findMaxPerturbation(stacks)
			perturbationLists[0].append(perturbations[0])
			perturbationLists[1].append(perturbations[1])
			perturbationLists[2].append(perturbations[2])
		perturbationMeshes[0].append(perturbationLists[0])
		perturbationMeshes[1].append(perturbationLists[1])
		perturbationMeshes[2].append(perturbationLists[2])

	ceresAngleLocation = simulationDataDirectory\
		+ prefix + "Sampled_Angles_from_Ceres.csv"
	ceresAngles = tools.transpose(
		tools.floatParseCSV(ceresAngleLocation))[0]
	for i in range(len(ceresAngles)):
		ceresAngles[i] *= 180 / np.pi
		# if ceresAngles[i] < 0:
		# 	ceresAngles[i] += 360

	for i in range(1, len(ceresAngles)):
		# print(ceresAngles[i], "  ", end = "")
		while abs(ceresAngles[i] - ceresAngles[i-1]) > 180:
			ceresAngles[i] += 360
		# print(ceresAngles[i])

	constellationAngleLocation = simulationDataDirectory\
		+ prefix + "Sampled_Constellation_Angles.csv"
	constellationAnglesFromFile = tools.transpose(
		tools.floatParseCSV(constellationAngleLocation))[0]
	names = ["Arm_1", "Arm_2", "Arm_3"]
	for i in range(3):
		graphMesh(ceresAngles, constellationAnglesFromFile,
			perturbationMeshes[i],
			prefix + names[i])

