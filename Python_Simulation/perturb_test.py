import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import my_package.util as tools
from math import *

G = 6.67408e-11
mu = G * 1.989e30

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

startYear = 2029
startEpoch = LISAutils.yearToEpochMJD(startYear)

earth = LISAutils.orbit(radians(200.7), 0.01671, radians(0),
	radians(-11.261), radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58324.75, name = 'earth')

initialConstellationPhase = 0
# initialConstellationPhase = 2 * np.pi / 6
lisa = LISAutils.LISA(initialConstellationPhase,
	np.arctan2(earth.y(), earth.x()))

ceres = LISAutils.orbit(radians(352.23), 0.0755, radians(10.5935),
	radians(80.3099), radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = 8.958e20)

year = 365 * 24 * 3600
T = 6 * year
dt = 1e5
ceres.setTime(0)
times = np.arange(0, T, dt)
ceresPositions = []
for t in times:
	ceres.setTime(t)
	ceresPositions.append(ceres.pos())
ceres.setTime(0)


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
	ceres.setTime(time)
	otherPos = ceres.pos()
	displacement = pos - otherPos
	posSqr = sum(displacement**2)
	specMu = ceres.mass * G
	return -displacement * specMu / posSqr**(1.5)

def solarDifference(pos, deltaPos):
	#returns solarForce(pos + deltaPos) - solarForce(pos)
	#up to second order in epsilon
	posSqr = sum(pos**2)
	epsilon = 2 * sum(pos * deltaPos) / posSqr + sum(deltaPos**2) / posSqr
	temp = 3 / 2 * epsilon - 15 / 8 * epsilon**2
	temp2 = mu / posSqr**1.5
	output = pos * temp2 * temp
	output -= deltaPos * temp2 * (1 - temp)
	return output

def stateDerivative(state, time):
	outputState = State(state.vel, solarForce(state.pos),
		state.dVel, bodyAccel(state.pos + state.dPos, time)
		+ solarDifference(state.pos, state.dPos))
	return outputState

def buildUnit(satA, satB):
	units = []
	size = len(satA)
	for i in range(size):
		units.append((satB[i] - satA[i]) / np.sqrt(sum((satB[i] - satA[i])**2)))
	return np.array(units)

def arrayDot(vecListA, vecListB):
	output = []
	for i in range(len(vecListA)):
		output.append(sum(vecListA[i] * vecListB[i]))
	return np.array(output)

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
		output.append(((angleThere - angleHere)%(2 * np.pi) - np.pi) * 180 / np.pi)
	return np.array(output)

stacks = []
for i in range(3):
	initialState = State(lisa.sat(i).pos(), lisa.sat(i).vel())
	stacks.append(rk4(initialState, stateDerivative, 0, T, dt))

tempThing = guidingCenterTrajectory(stacks)[0]

for stateList in stacks:
	stateItem = stateList[0]
	print(stateItem.pos - tempThing)
	# print(testingThing)
	# print(state.pos)
	# print(testingThing - state.pos)

newTimes = []
for t in times:
	newTimes.append(t / year + startYear)
newTimes = np.array(newTimes)


relAngle = relativeAngleToCeres(stacks)
diff = rHatPhiHatGuidingCenterDifference(stacks)
armDif = armLengthDifferences(stacks)
angleDif = armAngleDifferneces(stacks)

plt.figure(1)
# plt.plot(newTimes, relAngle, label = 'angle to Ceres', linestyle=':')

# plt.figure(2)
plt.plot(newTimes, diff[0], label = 'r displacement')
# plt.plot(newTimes, diff[1], label = 'phi displacement')
# plt.plot(newTimes, diff[2], label = 'z displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title('Perturbation of LISA Guiding Center due to Ceres')
# plt.legend()
# plt.show()

plt.figure(2)
plt.plot(newTimes, armDif[0], label = 'arm 1 displacement')
plt.plot(newTimes, armDif[1], label = 'arm 2 displacement')
plt.plot(newTimes, armDif[2], label = 'arm 3 displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title('Perturbation of LISA Arm Lengths due to Ceres\nwith 2 pi / 6 constellation phase angle')
plt.legend()
# plt.show()

plt.figure(3)
plt.plot(newTimes, angleDif[0] * 180 / np.pi * 3600, label = 'angle 1 displacement')
plt.plot(newTimes, angleDif[1] * 180 / np.pi * 3600, label = 'angle 2 displacement')
plt.plot(newTimes, angleDif[2] * 180 / np.pi * 3600, label = 'angle 3 displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('angle (arcseconds)')
plt.title('Perturbation of LISA Arm Angles due to Ceres')
plt.legend()
plt.show()

# sat1 = lisa.sat(0)

# state = []
# state.append(sat1.pos())
# state.append(sat1.vel())
# state.append(np.zeros((3)))
# state.append(np.zeros((3)))
# state = np.array(state)

# stateList = tools.transpose(rk4(state, stateDerivative, 0, T, dt))
# pos = stateList[0]
# pos = tools.transpose(pos)
# xPos = pos[0]
# yPos = pos[1]
# diff = stateList[2]
# diff = tools.transpose(diff)
# xList = diff[0]
# yList = diff[1]

# rHat = []
# phiHat = []
# for i in range(len(xList)):
# 	rHat.append((xPos[i] * xList[i] + yPos[i] * yList[i])
# 		/ np.sqrt(xPos[i] * xPos[i] + yPos[i] * yPos[i]))
# 	phiHat.append((-yPos[i] * xList[i] + xPos[i] * yList[i])
# 		/ np.sqrt(xPos[i] * xPos[i] + yPos[i] * yPos[i]))

# plt.plot(times / year, xList, label = 'x displacement')
# plt.plot(times / year, yList, label = 'y displacement')
# plt.plot(times / year, xPos / 2e8, label = 'x displacement')
# plt.plot(times / year, yPos / 2e8, label = 'y displacement')
# plt.plot(times / year, rHat, label = 'displacement in r direction')
# plt.plot(times / year, phiHat, label = 'displacement in phi direction')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('displacement (meters)')
# plt.title('Perturbation of LISA Satellite due to Ceres')
# # plt.plot(timeList, velFirstDif, label = 'vel first')
# # plt.plot(timeList, rk4Dif, label = 'RK4')
# plt.legend()
# plt.show()