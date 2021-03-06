import numpy as np
from Simulation_Core import *


def buildStacksAndTimes(fileLocation):
	data = tools.floatParseCSV(fileLocation)

	def buildState(satData, satPerturbedData):
		pos = np.array(satData[:3])
		vel = np.array(satData[3:])
		dPos = np.array(satPerturbedData[:3])
		dVel = np.array(satPerturbedData[3:])
		return State(pos, vel, dPos, dVel)

	stacks = [[],[],[]]
	times = []

	# print(len(data))

	for entry in data:
		#entry is a list of floats specifying time, lisa state, and delta lisa state
		times.append(entry[0])
		sat1Data = entry[1:7]
		sat2Data = entry[7:13]
		sat3Data = entry[13:19]
		sat1PerturbedData = entry[19:25]
		sat2PerturbedData = entry[25:31]
		sat3PerturbedData = entry[31:37]
		state1 = buildState(sat1Data, sat1PerturbedData)
		state2 = buildState(sat2Data, sat2PerturbedData)
		state3 = buildState(sat3Data, sat3PerturbedData)
		stacks[0].append(state1)
		stacks[1].append(state2)
		stacks[2].append(state3)

	return (stacks, times)

def getPos(stacks):
	output = []
	for item in stacks:
		temp = []
		for element in item:
			temp.append(element.pos)
		output.append(temp)
	return np.array(output)

def getD_Pos(stacks):
	output = []
	for item in stacks:
		temp = []
		for element in item:
			temp.append(element.dPos)
		output.append(temp)
	return np.array(output)

def getVel(stacks):
	output = []
	for item in stacks:
		temp = []
		for element in item:
			temp.append(element.vel)
		output.append(temp)
	return np.array(output)

def getD_Vel(stacks):
	output = []
	for item in stacks:
		temp = []
		for element in item:
			temp.append(element.dVel)
		output.append(temp)
	return np.array(output)

def armDifferencer(inputArray):
	return np.array([inputArray[2] - inputArray[1],
		inputArray[0] - inputArray[2],
		inputArray[1] - inputArray[0]])

def getArmVectors(stacks):
	return armDifferencer(getPos(stacks))

def getBaseArmLength(stacks):
	# armVectors = getArmVectors(stacks)
	armVectors = getArmVectors(stacks)
	output = []
	for arm in armVectors:
		armLengths = []
		for armVectorAtTime in arm:
			armLengths.append(np.sqrt(sum(armVectorAtTime * armVectorAtTime)))
		output.append(armLengths)
	return np.array(output)

def getArmLengths(stacks):
	# armVectors = getArmVectors(stacks)
	armVectors = getArmVectors(stacks) + armDifferencer(getD_Pos(stacks))
	output = []
	for arm in armVectors:
		armLengths = []
		for armVectorAtTime in arm:
			armLengths.append(np.sqrt(sum(armVectorAtTime * armVectorAtTime)))
		output.append(armLengths)
	return np.array(output)

def arrayNormalize(inputArray):
	return np.array([item / np.sqrt(sum(item * item)) for item in inputArray])

def getArmUnits(stacks):
	return np.array([arrayNormalize(item) for item in getArmVectors(stacks)])

def getD_ArmVectors(stacks):
	return armDifferencer(getD_Pos(stacks))

def arrayDot(inputArray1, inputArray2):
	return np.array([sum(item[0] * item[1]) for item in zip(inputArray1, inputArray2)])

def getD_ArmLengths(stacks):
	return np.array([arrayDot(item[0], item[1]) for item in zip(getD_ArmVectors(stacks), getArmUnits(stacks))])

def getArmAngles(stacks):
	armUnits = getArmUnits(stacks)
	return np.array([arrayDot(armUnits[1], -armUnits[2]),
		arrayDot(armUnits[2], -armUnits[0]),
		arrayDot(armUnits[0], -armUnits[1])])

def getD_ArmAngles(stacks):
	unperturbedAngles = getArmAngles(stacks)
	pos = getPos(stacks)
	dPos = getD_Pos(stacks)
	totalPos = pos + dPos
	newArmUnits = np.array([arrayNormalize(item) for item in armDifferencer(totalPos)])
	newAngles = np.array([arrayDot(newArmUnits[1], -newArmUnits[2]),
		arrayDot(newArmUnits[2], -newArmUnits[0]),
		arrayDot(newArmUnits[0], -newArmUnits[1])])
	return newAngles - unperturbedAngles

def getArmVel(stacks):
	vel = getVel(stacks)
	armVelVectors = armDifferencer(vel)
	return np.array([arrayDot(pair[0], pair[1]) for pair in zip(armVelVectors, getArmUnits(stacks))])

def getD_ArmVel(stacks):
	dVel = getD_Vel(stacks)
	armD_VelVectors = armDifferencer(dVel)
	return np.array([arrayDot(pair[0], pair[1]) for pair in zip(armD_VelVectors, getArmUnits(stacks))])

def getGuidingCenterPos(stacks):
	return np.array(sum(getPos(stacks))) / 3

def getGuidingCenterD_Pos(stacks):
	return np.array(sum(getD_Pos(stacks))) / 3

def get_rHatPhiHatGuidingCenterD_Pos(stacks):
	dPos = getGuidingCenterD_Pos(stacks)
	pos = getGuidingCenterPos(stacks)
	rUnit = arrayNormalize(pos)
	phiHat = []
	rHat = []
	for unit in rUnit:
		phiHat.append([-unit[1], unit[0], 0])
		rHat.append([unit[0], unit[1], 0])
	phiHat = np.array(phiHat)
	rHat = np.array(rHat)
	return np.array([arrayDot(rHat, dPos),
		arrayDot(phiHat, dPos)])

def getD_DistToGuide(stacks):
	guidePos = getGuidingCenterPos(stacks)
	guideD_Pos = getGuidingCenterD_Pos(stacks)
	output = []
	for item in zip(getPos(stacks), getD_Pos(stacks)):
		guideToPosHat = arrayNormalize(item[0] - guidePos)
		output.append(arrayDot(guideToPosHat, item[1] - guideD_Pos))
	return np.array(output)

def getAccels(stacks, times):
	output = []
	for satPos in getPos(stacks):
		temp = []
		for timePoint in zip(satPos, times):
			temp.append(bodyAccel(timePoint[0], timePoint[1]))
		output.append(temp)
	return np.array(output)

def getArmAccels(stacks, times, inCeres):
	ceresPositions = inCeres
	armUnits = getArmUnits(stacks)
	accels = armDifferencer(getAccels(stacks, times))
	return np.array([arrayDot(item[0], item[1]) for item in zip(armUnits, accels)])

def distToCeres(stacks, inCeres):
	guidePos = getGuidingCenterPos(stacks)
	output = []
	for timePoint in zip(inCeres, getGuidingCenterPos(stacks)):
		displacement = timePoint[0] - timePoint[1]
		output.append(np.sqrt(sum(displacement**2)))
	return np.array(output)
# def distToCeres(stacks, times, dt, inCeres):
# 	guidePos = getGuidingCenterPos(stacks)
# 	output = []
# 	for timePoint in zip(times, getGuidingCenterPos(stacks)):
# 		displacement = tools.linInterp(inCeres, 0, dt, timePoint[0]) - timePoint[1]
# 		output.append(np.sqrt(sum(displacement**2)))
# 	return np.array(output)

def getPerturbativeAccelValues(stacks):
	accels = [solarDifference(item[0], item[1]) for item in zip(getPos(stacks), getD_Pos(stacks))]
	return np.array([arrayDot(item[0], item[1]) for item in zip(armDifferencer(accels), getArmUnits(stacks))])

def getArmAccelDotVel(stacks):
	accels = getPerturbativeAccelValues(stacks) + getArmAccels(stacks)
	return np.array([arrayDot(item[0], item[1]) for item in zip(accels, getD_ArmVel(stacks))])

def normAngleFromArmsToCeres(stacks, times, dt, inCeres):
	armUnits = getArmUnits(stacks)
	dispToCeres = []
	for timePoint in zip(times, getGuidingCenterPos(stacks)):
		dispToCeres.append(tools.linInterp(inCeres, 0, dt, timePoint[0]) - timePoint[1])
	unitToCeres = arrayNormalize(dispToCeres)
	return np.array([np.arccos(arrayDot(armUnit, unitToCeres)) - np.pi / 2 for armUnit in armUnits])

def constellationArea(stacks):
	arms = getArmLengths(stacks)
	areas = []
	for armState in tools.transpose(arms):
		s = sum(armState)/2
		areas.append(np.sqrt(s*(s-armState[0])*(s-armState[1])*(s-armState[2])))
	return np.array(areas)

def getPerburbationOfArea(stacks):
	base = getBaseArmLength(stacks)
	pert = getArmLengths(stacks)
	dAreas = []
	for armState in zip(tools.transpose(base), tools.transpose(pert)):
		bs = sum(armState[0]) / 2
		ps = sum(armState[1]) / 2
		baseArea = np.sqrt(bs*(bs-armState[0][0])*(bs-armState[0][1])*(bs-armState[0][2]))
		pertArea = np.sqrt(ps*(ps-armState[1][0])*(ps-armState[1][1])*(ps-armState[1][2]))
		dAreas.append(pertArea - baseArea)
	return np.array(dAreas)

# def relativeAngleToCeres(inputStacks):
# 	pos = getGuidingCenterPos(inputStacks)
# 	output = []
# 	for i in range(len(pos)):
# 		angleHere = np.arctan2(pos[i][1], pos[i][0])
# 		angleThere = np.arctan2(ceresPositions[i][1], ceresPositions[i][0])
# 		output.append(((angleThere - angleHere + np.pi)%(2 * np.pi) - np.pi) * 180 / np.pi)
# 	return np.array(output)
