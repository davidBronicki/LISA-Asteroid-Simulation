import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# import matplotlib.colors as cl
import my_package.util as tools
from math import *

G = 6.67408e-11
mu = G * 1.989e30
year = 365 * 24 * 3600
dt = 0
ceresMass = 8.958e20
times = []
newTimes = []
ceresPositions = []

def rk4(initialState, stateDerivative, t0, t1, dt):
	timeList = np.arange(t0, t1 - dt, dt)
	state = initialState.copy()
	stateList = [state.copy()]
	for t in timeList:
		ds1 = dt * stateDerivative(state, t)
		ds2 = dt * stateDerivative(state + ds1 / 2, t + dt / 2)
		ds3 = dt * stateDerivative(state + ds2 / 2, t + dt / 2)
		ds4 = dt * stateDerivative(state + ds3, t + dt)
		state += (ds1 + 2 * ds2 + 2 * ds3 + ds4) / 6
		stateList.append(state.copy())
	return stateList

class State:
	def __init__(self, pos, vel, dPos = [], dVel = []):
		self.pos = pos
		self.vel = vel
		if len(dPos) == 0:
			self.dPos = np.zeros((len(pos)))
		else:
			self.dPos = dPos
		if len(dVel) == 0:
			self.dVel = np.zeros((len(pos)))
		else:
			self.dVel = dVel

	def copy(self):
		return State(self.pos.copy(), self.vel.copy(),
			self.dPos.copy(), self.dVel.copy())

	def __add__(self, other):
		return State(self.pos + other.pos,
			self.vel + other.vel,
			self.dPos + other.dPos,
			self.dVel + other.dVel)

	def __mul__(self, other):
		return State(self.pos * other,
			self.vel * other,
			self.dPos * other,
			self.dVel * other)

	def __rmul__(self, other):
		return self * other

	def __truediv__(self, other):
		return self * (1 / other)

	def __neg__(self):
		return self * -1


def solarForce(pos):
	posSqr = sum(pos * pos)
	return -pos * mu / posSqr**(1.5)

def bodyAccel(pos, time):
	otherPos = ceresPositions[int(time / dt + 0.5)]
	displacement = pos - otherPos
	posSqr = sum(displacement**2)
	specMu = ceresMass * G
	return -displacement * specMu / posSqr**(1.5)

def solarDifference(pos, deltaPos):
	#returns solarForce(pos + deltaPos) - solarForce(pos)
	#up to second order in epsilon (defined below)
	posSqr = sum(pos**2)
	epsilon = 2 * sum(pos * deltaPos) / posSqr + sum(deltaPos**2) / posSqr
	alpha = 3 / 2 * epsilon - 15 / 8 * epsilon**2
	coef = mu / posSqr**1.5
	output = pos * coef * alpha
	output -= deltaPos * coef * (1 - alpha)
	return output

def stateDerivative(state, time):
	outputState = State(state.vel, solarForce(state.pos),
		state.dVel, bodyAccel(state.pos + state.dPos, time)
		+ solarDifference(state.pos, state.dPos))
	return outputState







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

def getArmVectors(stacks):
	positions = getPos(stacks)
	return np.array([positions[2] - positions[1],
		positions[0] - positions[2],
		positions[1] - positions[0]])

def getArmLengths(stacks):
	armVectors = getArmVectors(stacks)
	output = []
	for arm in armVectors:
		armLengths = []
		for armVectorAtTime in arm:
			armLengths.append(np.sqrt(sum(armVectorAtTime * armVectorAtTime)))
		output.append(armLengths)
	return np.array(output)

def arrayNormalize(inputArray):
	output = []
	for vector in inputArray:
		output.append(vector / np.sqrt(sum(vector * vector)))
	return np.array(output)

def getArmUnits(stacks):
	armVectors = getArmVectors(stacks)
	output = []
	for armVector in armVectors:
		output.append(arrayNormalize(armVector))
	return np.array(output)

def getD_ArmVectors(stacks):
	dPositions = getD_Pos(stacks)
	return np.array([dPositions[2] - dPositions[1],
		dPositions[0] - dPositions[2],
		dPositions[1] - dPositions[0]])

def arrayDot(inputArray1, inputArray2):
	output = []
	for pair in zip(inputArray1, inputArray2):
		output.append(sum(pair[0] * pair[1]))
	return np.array(output)

def getD_ArmLengths(stacks):
	dArmVectors = getD_ArmVectors(stacks)
	armUnits = getArmUnits(stacks)
	output = []
	for pair in zip(dArmVectors, armUnits):
		output.append(arrayDot(pair[0], pair[1]))
	return np.array(output)

def getArmAngles(stacks):
	armUnits = getArmUnits(stacks)
	return np.array([arrayDot(armUnits[1], -armUnits[2]),
		arrayDot(armUnits[2], -armUnits[0]),
		arrayDot(armUnits[0], -armUnits[1])])

def getGuidingCenterPos(stacks):
	positions = getPos(stacks)
	return np.array(sum(positions)) / 3

def getGuidingCenterD_Pos(stacks):
	dPositions = getD_Pos(stacks)
	return np.array(sum(dPositions)) / 3

def getAccels(stacks):
	output = []
	for satPos in getPos(stacks):
		temp = []
		for timePoint in zip(satPos, times):
			temp.append(bodyAccel(timePoint[0], timePoint[1]))
		output.append(temp)
	return np.array(output)

def getArmAccels(stacks):
	armUnits = getArmUnits(stacks)

	positions = np.array([
		collectPos(stacks[0]),
		collectPos(stacks[1]),
		collectPos(stacks[2])])

	accels = [[],[],[]]

	for i in range(len(ceresPositions)):
		accels[0].append(bodyAccel(positions[0][i], times[i]))
		accels[1].append(bodyAccel(positions[1][i], times[i]))
		accels[2].append(bodyAccel(positions[2][i], times[i]))

	thing = [[],[],[]]

	thing[0] = arrayDot(np.array(accels[2]) - np.array(accels[1]), buildUnit(positions[1], positions[2]))
	thing[1] = arrayDot(np.array(accels[2]) - np.array(accels[0]), buildUnit(positions[0], positions[2]))
	thing[2] = arrayDot(np.array(accels[0]) - np.array(accels[1]), buildUnit(positions[1], positions[0]))

	return np.array(thing)









def vectorToUnit(inputVectorArray):
	output = []
	for item in inputVectorArray:
		newThing = item / sqrt(sum(item * item))
		output.append(newThing)
	return np.array(output)

def buildUnit(satA, satB):
	units = []
	size = len(satA)
	for i in range(size):
		units.append((satB[i] - satA[i]) / np.sqrt(sum((satB[i] - satA[i])**2)))
	return np.array(units)

def arrayNormalize(vecList):
	output = []
	for item in vecList:
		output.append(item / sum(item * item))
	return np.array(output)

# def arrayDot(vecListA, vecListB):
# 	output = []
# 	for i in range(len(vecListA)):
# 		output.append(sum(vecListA[i] * vecListB[i]))
# 	return np.array(output)

def collectPos(inputList):
	output = []
	for item in inputList:
		output.append(item.pos)
	return np.array(output)

def collect_dPos(inputList):
	output = []
	for item in inputList:
		output.append(item.dPos)
	return np.array(output)

def collectVel(inputList):
	output = []
	for item in inputList:
		output.append(item.vel)
	return np.array(output)

def collect_dVel(inputList):
	output = []
	for item in inputList:
		output.append(item.dVel)
	return np.array(output)

def guidingCenterTrajectory(inputStacks):
	result = np.zeros((len(inputStacks[0]), 3))
	for sat in inputStacks:
		for i in range(len(sat)):
			result[i] += sat[i].pos
	result /= 3
	return result

def guidingCenterDifference(inputStacks):
	result = np.zeros((len(inputStacks[0]), 3))
	for sat in inputStacks:
		for i in range(len(sat)):
			result[i] += sat[i].dPos
	result /= 3
	return result

def rHatPhiHatGuidingCenterDifference(inputStacks):
	pos = guidingCenterTrajectory(inputStacks)
	dif = guidingCenterDifference(inputStacks)
	rHat = []
	phiHat = []
	for spot in pos:
		tempVal = np.array(spot[0], spot[1])
		rHat.append(spot / sqrt(sum(spot * spot)))
		phiHat.append(np.array(-spot[1], spot[0]) / sqrt(sum(spot * spot)))
	return np.array([arrayDot(rHat, dif), arrayDot(phiHat, dif)])

def armLengthDifferences(inputStacks):
	sat1Pos = collectPos(inputStacks[0])
	sat2Pos = collectPos(inputStacks[1])
	sat3Pos = collectPos(inputStacks[2])

	unit12 = buildUnit(sat1Pos, sat2Pos)
	unit13 = buildUnit(sat1Pos, sat3Pos)
	unit23 = buildUnit(sat2Pos, sat3Pos)

	disp1 = collect_dPos(inputStacks[0])
	disp2 = collect_dPos(inputStacks[1])
	disp3 = collect_dPos(inputStacks[2])

	result = []
	result.append(arrayDot(disp3 - disp2, unit23))
	result.append(arrayDot(disp3 - disp1, unit13))
	result.append(arrayDot(disp2 - disp1, unit12))
	return np.array(result)

def armAngle(a, b, c):
	unit1 = buildUnit(b, a)
	unit2 = buildUnit(b, c)
	return np.arccos(arrayDot(unit1, unit2))

def armAngleDifferneces(inputStacks):
	sat1Pos = collectPos(inputStacks[0])
	sat2Pos = collectPos(inputStacks[1])
	sat3Pos = collectPos(inputStacks[2])

	disp1 = collect_dPos(inputStacks[0])
	disp2 = collect_dPos(inputStacks[1])
	disp3 = collect_dPos(inputStacks[2])

	pert1 = sat1Pos + disp1
	pert2 = sat2Pos + disp2
	pert3 = sat3Pos + disp3

	output = []
	output.append(armAngle(pert2, pert1, pert3) - armAngle(sat2Pos, sat1Pos, sat3Pos))
	output.append(armAngle(pert1, pert2, pert3) - armAngle(sat1Pos, sat2Pos, sat3Pos))
	output.append(armAngle(pert1, pert3, pert2) - armAngle(sat1Pos, sat3Pos, sat2Pos))
	return np.array(output)

def relativeAngleToCeres(inputStacks):
	pos = guidingCenterTrajectory(inputStacks)
	output = []
	for i in range(len(pos)):
		angleHere = np.arctan2(pos[i][1], pos[i][0])
		angleThere = np.arctan2(ceresPositions[i][1], ceresPositions[i][0])
		output.append(((angleThere - angleHere + np.pi)%(2 * np.pi) - np.pi) * 180 / np.pi)
	return np.array(output)

def normAngleFromArmsToCeres(inputStacks):
	sat1Pos = collectPos(inputStacks[0])
	sat2Pos = collectPos(inputStacks[1])
	sat3Pos = collectPos(inputStacks[2])

	unit12 = buildUnit(sat1Pos, sat2Pos)
	unit13 = buildUnit(sat1Pos, sat3Pos)
	unit23 = buildUnit(sat2Pos, sat3Pos)

	dispToCeres = ceresPositions - sat1Pos
	unitToCeres = vectorToUnit(dispToCeres)

	arm1DotVecToCeres = arrayDot(unitToCeres, unit23)
	arm2DotVecToCeres = arrayDot(unitToCeres, unit13)
	arm3DotVecToCeres = arrayDot(unitToCeres, unit12)

	return np.array([np.arccos(arm1DotVecToCeres) - np.pi / 2,
		np.arccos(arm2DotVecToCeres) - np.pi / 2,
		np.arccos(arm3DotVecToCeres) - np.pi / 2])

def diffRelPosToGuidingCenter(inputStacks):
	output = []

	guidingCenter = guidingCenterTrajectory(inputStacks)
	guidingCenterDiff = guidingCenterDifference(inputStacks)

	sat1Pos = collectPos(inputStacks[0])
	guideTo1Vec = sat1Pos - guidingCenter
	guideTo1Hat = vectorToUnit(guideTo1Vec)
	dispSat1 = collect_dPos(inputStacks[0])
	guideTo1Diff = dispSat1 - guidingCenterDiff
	deltaDistToGuide1 = arrayDot(guideTo1Diff, guideTo1Hat)
	output.append(deltaDistToGuide1)

	sat2Pos = collectPos(inputStacks[1])
	guideTo2Vec = sat2Pos - guidingCenter
	guideTo2Hat = vectorToUnit(guideTo2Vec)
	dispSat2 = collect_dPos(inputStacks[1])
	guideTo2Diff = dispSat2 - guidingCenterDiff
	deltaDistToGuide2 = arrayDot(guideTo2Diff, guideTo2Hat)
	output.append(deltaDistToGuide2)

	sat3Pos = collectPos(inputStacks[2])
	guideTo3Vec = sat3Pos - guidingCenter
	guideTo3Hat = vectorToUnit(guideTo3Vec)
	dispSat3 = collect_dPos(inputStacks[2])
	guideTo3Diff = dispSat3 - guidingCenterDiff
	deltaDistToGuide3 = arrayDot(guideTo3Diff, guideTo3Hat)
	output.append(deltaDistToGuide3)

	return np.array(output)

def getPerturbativeAccelValues(stacks):
	pos1 = collectPos(stacks[0])
	pos2 = collectPos(stacks[1])
	pos3 = collectPos(stacks[2])

	dPos1 = collect_dPos(stacks[0])
	dPos2 = collect_dPos(stacks[1])
	dPos3 = collect_dPos(stacks[2])

	accel1 = []
	accel2 = []
	accel3 = []

	for i in range(len(pos1)):
		accel1.append(solarDifference(pos1[i], dPos1[i]))
		accel2.append(solarDifference(pos2[i], dPos2[i]))
		accel3.append(solarDifference(pos3[i], dPos3[i]))

	output = [[],[],[]]

	output[0] = arrayDot(np.array(accel3) - np.array(accel2), buildUnit(pos2, pos3))
	output[1] = arrayDot(np.array(accel3) - np.array(accel1), buildUnit(pos1, pos3))
	output[2] = arrayDot(np.array(accel2) - np.array(accel1), buildUnit(pos1, pos2))

	return np.array(output)

def getArmVelDiff(stacks):
	positions = np.array([
		collectPos(stacks[0]),
		collectPos(stacks[1]),
		collectPos(stacks[2])])

	velocities = np.array([
		collect_dVel(stacks[0]),
		collect_dVel(stacks[1]),
		collect_dVel(stacks[2])])

	arm_dVels = np.array([
		arrayDot(velocities[0] - velocities[1], buildUnit(positions[1], positions[0])),
		arrayDot(velocities[1] - velocities[2], buildUnit(positions[2], positions[1])),
		arrayDot(velocities[2] - velocities[0], buildUnit(positions[0], positions[2]))])

	return arm_dVels

	# for i in range(len(armLengths)):
	# 	for j in range(len(armLengths[0])):
	# 		armLengths[i][j] = np.sqrt(armLengths[i][j])

	# return armLengths

def buildArmAccelDotVel(stacks):
	accels = getPerturbativeAccelValues(stacks) + getArmAccels(stacks)
	vels = getArmVelDiff(stacks)

	output = [[],[],[]]
	for i in range(len(vels[0])):
		output[0].append(accels[0][i] * vels[0][i])
		output[1].append(accels[1][i] * vels[1][i])
		output[2].append(accels[2][i] * vels[2][i])

	return np.array(output)

def performSimulation(initialConstellationPhase, startYear, totalTime, inputDT):
	global times, ceresPositions, newTimes, dt
	times = []
	ceresPositions = []
	newTimes = []
	dt = inputDT

	startEpoch = LISAutils.yearToEpochMJD(startYear)

	earth = LISAutils.orbit(radians(200.7), 0.01671, radians(0),
		radians(-11.261), radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58324.75, name = 'earth')
	lisa = LISAutils.LISA(initialConstellationPhase,
		np.arctan2(earth.y(), earth.x()))
		# np.pi * 1/6 + np.pi * 2 / 3)
	ceres = LISAutils.orbit(radians(352.23), 0.0755, radians(10.5935),
		radians(80.3099), radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = ceresMass)

	# print(earth.pos())
	# print(lisa.pos())

	times = np.arange(0, totalTime, dt)
	for t in times:
		ceres.setTime(t)
		ceresPositions.append(ceres.pos())
		newTimes.append(t / year + startYear)
		# print(newTimes[-1])
	newTimes = np.array(newTimes)
	ceresPositions = np.array(ceresPositions)
	ceres.setTime(0)

	stacks = []
	for i in range(3):
		initialState = State(lisa.sat(i).pos(), lisa.sat(i).vel())
		stacks.append(rk4(initialState, stateDerivative, 0, totalTime, dt))

	# print(len(stacks[0]))
	# print(len(ceresPositions))
	# print(len(times))

	return stacks


stacks = performSimulation(np.pi / 6 * 0, 2029, 6*year, 1e5)

# relAngle = relativeAngleToCeres(stacks)
# normAngleToCeres = normAngleFromArmsToCeres(stacks)
diff = getGuidingCenterD_Pos(stacks)
armDif = getD_ArmLengths(stacks)
# armDif = armLengthDifferences(stacks)
# deltaDistToGuide = getGuidingCenterD_Pos(stacks)
armAccels = getArmAccels(stacks)
# twoBodyPerturbAccels = getPerturbativeAccelValues(stacks)
armLengths = getArmLengths(stacks)
# angleDif = armAngleDifferneces(stacks)
# armAccelDotVel = buildArmAccelDotVel(stacks)



##############build mp4 of arm perturbation with varying constellation angles###############

# fig, ax = plt.subplots()

# x = np.arange(0, 2*np.pi, 0.01)
# plt.ylim([-18, 18])
# plt.xlabel('time (year)')
# plt.ylabel('displacement (meters)')
# plt.title('Perturbation of LISA Arm 1 due to Ceres')
# plt.grid()
# line, = ax.plot(newTimes, armLengthDifferences(buildStacks(0))[0], label = 'phi = 0')

# frameCount = 180

# stacksList = []
# for i in range(frameCount):
# 	stacksList.append(buildStacks(i * 2 * np.pi / frameCount))
# 	print('initializing stacks: ', i * 2 * np.pi / frameCount)

# def animate(frame):
# 	line.set_ydata(armLengthDifferences(stacksList[frame])[0])
# 	line.set_label('phi = ' + str(round(frame * 2 * np.pi / frameCount, 2)))
# 	plt.legend()
# 	print('building frames: ', frame * 2 * np.pi / frameCount)
# 	return line,


# animation = ani.FuncAnimation(
#     fig, animate, interval=100, save_count=frameCount)

# animation.save("thing.mp4")

#########displays offset from center for pi/3 validation#######

# tempThing = guidingCenterTrajectory(stacks)[0]

# for stateList in stacks:
# 	stateItem = stateList[0]
# 	print(stateItem.pos - tempThing)
# 	# print(testingThing)
# 	# print(state.pos)
# 	# print(testingThing - state.pos)

########################Comparison to Analytic Results########################


# def rcpX(alpha, beta, e):
# 	return LISAutils.astroUnit * (np.cos(alpha)
# 		+ .5 * e * (np.cos(2 * alpha - beta) - 3 * np.cos(beta))
# 		+ .125 * e**2 * (3 * np.cos(3 * alpha - 2 * beta)
# 			- 10 * np.cos(alpha) - 5 * np.cos(alpha - 2 * beta)))

# def rcpY(alpha, beta, e):
# 	return LISAutils.astroUnit * (np.sin(alpha)
# 		+ .5 * e * (np.sin(2 * alpha - beta) - 3 * np.sin(beta))
# 		+ .125 * e**2 * (3 * np.sin(3 * alpha - 2 * beta)
# 			- 10 * np.sin(alpha) + 5 * np.sin(alpha - 2 * beta)))

# def rcpZ(alpha, beta, e):
# 	return LISAutils.astroUnit * np.sqrt(3) * e * (-cos(alpha - beta)
# 		+ e * (cos(alpha - beta)**2 + 2 * sin(alpha - beta)**2))

# def rcpPos(t, satNumber, initialConstellationPhase, initialLongitude):
# 	omega = sqrt(mu/LISAutils.astroUnit**3)

# 	alpha = omega*t + initialLongitude
# 	beta = 2 * np.pi * satNumber / 3 + initialConstellationPhase
# 	temp = 2.5e9 / (2 * LISAutils.astroUnit)
# 	e = np.sqrt(1 + (2 / np.sqrt(3)) * temp + (4 / 3) * temp ** 2) - 1

# 	return np.array([rcpX(alpha, beta, e),
# 		rcpY(alpha, beta, e), rcpZ(alpha, beta, e)])

# def rcpArmLengths(t, initialConstellationPhase, initialLongitude):
# 	omega = sqrt(mu/LISAutils.astroUnit**3)
# 	alpha = omega*t + initialLongitude
# 	temp = 2.5e9 / (2 * LISAutils.astroUnit)
# 	e = np.sqrt(1 + (2 / np.sqrt(3)) * temp + (4 / 3) * temp ** 2) - 1
# 	return 2.5e9 * np.array([
# 		1 - e/32 * (15 * cos(alpha - initialConstellationPhase)
# 			+ cos(3 * (alpha - initialConstellationPhase))),
# 		1 + e/32 * (-15 * sin(alpha - initialConstellationPhase - np.pi / 6)
# 			- cos(3 * (alpha - initialConstellationPhase))),
# 		1 + e/32 * (15 * sin(alpha - initialConstellationPhase + np.pi / 6)
# 			- cos(3 * (alpha - initialConstellationPhase)))])


# intermediateUnitlessVal = 2.5e9 / (2.0 * LISAutils.astroUnit)
# lisaEccentricity = np.sqrt(1.0 + (2.0 / sqrt(3.0)) * intermediateUnitlessVal
# 	+ (4.0 / 3.0) * intermediateUnitlessVal**2.0) - 1.0
# lisaMeanAnglularVelocity = sqrt(mu / LISAutils.astroUnit**3)

# def alphaFunct(t, initialLongitude, initialConstellationPhase):
#     return t * lisaMeanAnglularVelocity + initialLongitude - initialConstellationPhase

# def betaFunct(satNumber, initialConstellationPhase, initialLongitude):
#     return 2.0 * np.pi * satNumber / 3.0 + initialConstellationPhase - initialLongitude

# def psiFunct(alpha,beta):
#     return (alpha - beta - lisaEccentricity * sin(alpha - beta)
#     	+ (lisaEccentricity**2.0) * cos(alpha - beta) * sin(alpha - beta))

# def x(psi, satNumber):
#     x = ((LISAutils.astroUnit * (cos(psi) + lisaEccentricity) * (1.0 + intermediateUnitlessVal / sqrt(3.0))
# 	    	/ (1.0 + lisaEccentricity)) * cos((2.0 * np.pi / 3.0) * (satNumber - 1)) 
# 	    - (LISAutils.astroUnit * sqrt(1.0 - lisaEccentricity**2.0) * sin(psi)) * sin((2.0 * np.pi / 3.0) * (satNumber - 1)))
#     return x

# def y(psi, satNumber):
#     y = ((LISAutils.astroUnit * (cos(psi) + lisaEccentricity) * (1.0 + intermediateUnitlessVal / sqrt(3.0))
# 	    	/ (1.0 + lisaEccentricity)) * sin((2.0 * np.pi / 3.0) * (satNumber - 1))
# 	    + (LISAutils.astroUnit * sqrt(1.0 - lisaEccentricity**2.0) * sin(psi)) * cos((2.0 * np.pi / 3.0) * (satNumber - 1)))
#     return y

# def z(psi):
#     z = LISAutils.astroUnit * (cos(psi) + lisaEccentricity) * (intermediateUnitlessVal / (1.0 + lisaEccentricity))
#     return z

# def dnkvPos(t, satNumber, initialConstellationPhase, initialLongitude):
# 	psi = psiFunct(alphaFunct(t, initialLongitude, initialConstellationPhase),
# 		betaFunct(satNumber, initialConstellationPhase, initialLongitude))
# 	return np.array([
# 		x(psi, satNumber), y(psi, satNumber), z(psi)])

# rcpArms = []
# dnkvArms = [[],[],[]]
# for t in np.arange(0, 6*year, 1e4):
# 	# earth = LISAutils.orbit(radians(200.7), 0.01671, radians(0),
# 	# 	radians(-11.261), radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
# 	# 	simTime = LISAutils.yearToEpochMJD(2029), paramTime = 58324.75, name = 'earth')
# 	pos1 = dnkvPos(t, 1, 0, 0)#np.arctan2(earth.y(), earth.x()) - np.pi / 9)
# 	pos2 = dnkvPos(t, 2, 0, 0)#np.arctan2(earth.y(), earth.x()) - np.pi / 9)
# 	pos3 = dnkvPos(t, 3, 0, 0)#np.arctan2(earth.y(), earth.x()) - np.pi / 9)
# 	dnkvArms[0].append(np.sqrt(sum((pos1 - pos2) * (pos1 - pos2))))
# 	dnkvArms[1].append(np.sqrt(sum((pos1 - pos3) * (pos1 - pos3))))
# 	dnkvArms[2].append(np.sqrt(sum((pos3 - pos2) * (pos3 - pos2))))
# 	rcpArms.append(rcpArmLengths(t, 0, 0))#np.arctan2(earth.y(), earth.x())))

# rcpArms = tools.transpose(rcpArms)
# dnkvArms = np.array(dnkvArms)

# plt.figure(1)
# plt.plot(newTimes, rcpArms[0])
# plt.plot(newTimes, rcpArms[1])
# plt.plot(newTimes, rcpArms[2])
# plt.figure(2)
# plt.plot(newTimes, dnkvArms[0])
# plt.plot(newTimes, dnkvArms[1])
# plt.plot(newTimes, dnkvArms[2])
# plt.figure(3)
# plt.plot(newTimes, armLengths[0])
# plt.plot(newTimes, armLengths[1])
# plt.plot(newTimes, armLengths[2])
# plt.figure(4)
# plt.plot(newTimes, armLengths[0] - dnkvArms[0])
# plt.plot(newTimes, armLengths[1] - dnkvArms[1])
# plt.plot(newTimes, armLengths[2] - dnkvArms[2])
# plt.title('Difference from Integrator to DNKV Analytic Result\n1000 Second Step Size')
# plt.xlabel('Time (year)')
# plt.ylabel('Displacement (meters)')
# plt.grid()
# plt.show()



################energy argument plots#####################

# plt.figure(1)
# plt.plot(newTimes, twoBodyPerturbAccels[0])
# plt.plot(newTimes, twoBodyPerturbAccels[1])
# plt.plot(newTimes, twoBodyPerturbAccels[2])
# plt.grid()
# plt.title('perturbative acceleration')
# plt.figure(2)
# plt.plot(newTimes, armAccels[0])
# plt.plot(newTimes, armAccels[1])
# plt.plot(newTimes, armAccels[2])
# plt.grid()
# plt.title('ceres acceleration')
# plt.figure(3)
# plt.plot(newTimes, twoBodyPerturbAccels[0] + armAccels[0])
# plt.plot(newTimes, twoBodyPerturbAccels[1] + armAccels[1])
# plt.plot(newTimes, twoBodyPerturbAccels[2] + armAccels[2])
# plt.grid()
# plt.title('total acceleration')
plt.figure(4)
plt.plot(newTimes, armDif[0])
plt.plot(newTimes, armDif[1])
plt.plot(newTimes, armDif[2])
plt.grid()
plt.title('arm perturbation')
# plt.figure(5)
# plt.plot(newTimes, armAccelDotVel[0])
# plt.plot(newTimes, armAccelDotVel[1])
# plt.plot(newTimes, armAccelDotVel[2])
# plt.grid()
# plt.title('arm acceleration dot velocity')

# tempThing1 = 0
# finalThing1 = []
# tempThing2 = 0
# finalThing2 = []
# tempThing3 = 0
# finalThing3 = []
# for i in range(len(armAccelDotVel[0])):
# 	tempThing1 += armAccelDotVel[0][i] * 1e5
# 	finalThing1.append(tempThing1)
# 	tempThing2 += armAccelDotVel[1][i] * 1e5
# 	finalThing2.append(tempThing2)
# 	tempThing3 += armAccelDotVel[2][i] * 1e5
# 	finalThing3.append(tempThing3)

# plt.figure(6)
# plt.plot(newTimes, finalThing1)
# plt.plot(newTimes, finalThing2)
# plt.plot(newTimes, finalThing3)
# plt.grid()
# plt.title('energy')

plt.show()

##########multi-plot of different initial constellation phases###################


# def buildColor(angle):
# 	if angle < 40:
# 		theta = angle * (9 / 4) * np.pi / 180
# 		return (np.cos(theta), np.sin(theta), 0, 1)
# 	elif angle < 80:
# 		theta = (angle - 40) * (9 / 4) * np.pi / 180
# 		return (0, np.cos(theta), np.sin(theta), 1)
# 	else:
# 		theta = (angle - 80) * (9 / 4) * np.pi / 180
# 		return (np.sin(theta), 0, np.cos(theta), 1)

# # plt.rc('text', usetex=True)
# # plt.rc('font', family='serif')
# # for angle in range(0, 185, 20):
# for angle in range(0, 60, 5):
# 	theta = angle * np.pi / 180
# 	stacks = performSimulation(theta, 2029, 6 * year, 1e6)
# 	armAccels = getArmAccels(stacks)

# 	# plt.plot(newTimes, armAccels[0] * 1e15, color = buildColor(angle * 120 / 180))
# 	plt.figure(1)
# 	plt.plot(newTimes, armAccels[0] * 1e15, color = buildColor(2 * angle * 2 /3))
# 	plt.figure(2)
# 	plt.plot(newTimes, armAccels[1] * 1e15, color = buildColor(2 * angle * 2 /3))
# 	plt.figure(3)
# 	plt.plot(newTimes, armAccels[2] * 1e15, color = buildColor(2 * angle * 2 /3))
# 	# print(angle)
# 	# plt.figure(2)
# 	# normAngleToCeres = normAngleFromArmsToCeres(stacks)
# 	# plt.plot(newTimes, normAngleToCeres[0]*180/np.pi,# label = 'arm 1 normal angle to Ceres',
# 	# 	color = buildColor(2 * angle * 2 / 3))
# # plt.legend()
# plt.figure(1)
# plt.grid()
# plt.title('Arm 1 Acceleration due to Ceres')
# plt.xlabel('time (year)')
# plt.ylabel('acceleration (fm/s^2)')
# plt.ylim([-10, 10])
# plt.figure(2)
# plt.grid()
# plt.title('Arm 2 Acceleration due to Ceres')
# plt.xlabel('time (year)')
# plt.ylabel('acceleration (fm/s^2)')
# plt.ylim([-10, 10])
# plt.figure(3)
# plt.grid()
# plt.title('Arm 3 Acceleration due to Ceres')
# plt.xlabel('time (year)')
# plt.ylabel('acceleration (fm/s^2)')
# plt.ylim([-10, 10])
# # plt.show()


#####################Standard plotting routines################################


# plt.figure(1)
# plt.plot(newTimes, deltaDistToGuide[0], label = "Satellite 1")
# plt.plot(newTimes, deltaDistToGuide[1], label = "Satellite 2")
# plt.plot(newTimes, deltaDistToGuide[2], label = "Satellite 3")
# plt.grid()
# plt.title('Perturbation of the Distance of LISA Satellites\nto the Guiding Center (Constellation Angle 120)')
# plt.xlabel('time (year)')
# plt.ylabel('displacement (meters)')
# # plt.ylim([-18, 18])
# plt.legend()


plt.figure(1)
# plt.plot(newTimes, relAngle, label = 'angle to Ceres', linestyle=':')
plt.plot(newTimes, diff, label = 'r displacement')
# plt.plot(newTimes, diff[1], label = 'phi displacement')
# plt.plot(newTimes, diff[2], label = 'z displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title('Perturbation of LISA Guiding Center due to Ceres\nin r-hat Direction')
# plt.ylim([-180, 180])
# plt.legend()

plt.figure(2)
plt.plot(newTimes, armDif[0], label = 'arm 1 displacement')
plt.plot(newTimes, armDif[1], label = 'arm 2 displacement')
plt.plot(newTimes, armDif[2], label = 'arm 3 displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title('Perturbation of LISA Arm Lengths due to Ceres')#\nwith 0 constellation phase angle')
plt.legend()
# plt.ylim([-30, 30])

# plt.figure(3)
# # plt.plot(newTimes, normAngleToCeres[0]*180/np.pi, label = 'arm 1 normal angle to Ceres')
# # plt.plot(newTimes, normAngleToCeres[1]*180/np.pi, label = 'arm 2 normal angle to Ceres')
# # plt.plot(newTimes, normAngleToCeres[2]*180/np.pi, label = 'arm 3 normal angle to Ceres')
# plt.plot(newTimes, relAngle, color = 'black')# label = 'half angle from guiding center to Ceres', color = 'black')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('angle (degrees)')
# plt.title('Relative Angles Between Normal to LISA Arm 1 and Ceres. Black line is angle from LISA Constellation to Ceres\nColor Mapped According to Constellation Angle on x-axis: red = 0 degrees, blue = 60 degrees')#\nwith 0 constellation phase angle')
# # plt.legend()
# plt.ylim([-90, 90])

plt.figure(4)
plt.plot(newTimes, 10**15*armAccels[0], label = 'arm 1 acceleration')
plt.plot(newTimes, 10**15*armAccels[1], label = 'arm 2 acceleration')
plt.plot(newTimes, 10**15*armAccels[2], label = 'arm 3 acceleration')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('acceleration (pN/kg)')
plt.title('Perturbative Acceleration due to Ceres')
plt.legend()
plt.ylim([-12, 12])

# # plt.figure(5)
# # plt.plot(newTimes, angleDif[0] * 180 / np.pi * 3600, label = 'angle 1 displacement')
# # plt.plot(newTimes, angleDif[1] * 180 / np.pi * 3600, label = 'angle 2 displacement')
# # plt.plot(newTimes, angleDif[2] * 180 / np.pi * 3600, label = 'angle 3 displacement')
# # plt.grid()
# # plt.xlabel('time (years)')
# # plt.ylabel('angle (arcseconds)')
# # plt.title('Perturbation of LISA Arm Angles due to Ceres')
# # plt.legend()
# # plt.ylim([-0.0015, 0.0015])
plt.show()