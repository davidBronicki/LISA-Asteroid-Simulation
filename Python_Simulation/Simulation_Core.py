import numpy as np
import LISAutils
import math
import my_package.util as tools
import time as timePackage

G = 6.67408e-11
mu = G * 1.989e30
year = 365 * 24 * 3600
ceresMass = 8.958e20
ceresPositions = []
dt = 0
epsilon = 1e-10

def rk4Step(state, stateDerivative, t, dt):
	ds1 = dt * stateDerivative(state, t)
	ds2 = dt * stateDerivative(state + ds1 / 2, t + dt / 2)
	ds3 = dt * stateDerivative(state + ds2 / 2, t + dt / 2)
	ds4 = dt * stateDerivative(state + ds3, t + dt)
	return state + (ds1 + 2 * ds2 + 2 * ds3 + ds4) / 6

def rk4(initialState, stateDerivative, numericalErrorFunction, t0, t1, init_dt):
	dt = init_dt
	t = t0
	timeList = [t]
	stateList = [initialState.copy()]
	while (t < t1):
		timeList.append(t + dt)
		timeList.append(t + 2 * dt)
		testState = rk4Step(stateList[-1], stateDerivative, t, 2 * dt)
		stateList.append(rk4Step(stateList[-1], stateDerivative, t, dt))
		stateList.append(rk4Step(stateList[-1], stateDerivative, t + dt, dt))
		t += 2 * dt
		deltaState = testState - stateList[-1]
		ds = numericalErrorFunction(stateList[-1], deltaState, dt)
		if (ds < 1):
			dt *= ds**(1/4)
		else:
			dt *= ds**(1/5)
		# print(dt)
	return timeList, stateList

# def rk4(initialState, stateDerivative, t0, t1, dt):
# 	timeList = np.arange(t0, t1 - dt, dt)
# 	state = initialState.copy()
# 	stateList = [state.copy()]
# 	for t in timeList:
# 		ds1 = dt * stateDerivative(state, t)
# 		ds2 = dt * stateDerivative(state + ds1 / 2, t + dt / 2)
# 		ds3 = dt * stateDerivative(state + ds2 / 2, t + dt / 2)
# 		ds4 = dt * stateDerivative(state + ds3, t + dt)
# 		state += (ds1 + 2 * ds2 + 2 * ds3 + ds4) / 6
# 		stateList.append(state.copy())
# 	return stateList

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

	def toArray(self):
		return np.array([self.pos, self.vel, self.dPos, self.dVel])

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

	def __sub__(self, other):
		return self + (-other)

class LISA_State:
	def __init__(self, initialConstellationPhase, initialLongitude, sat3 = None):
		if sat3 != None:
			self.sat1 = initialConstellationPhase.copy()
			self.sat2 = initialLongitude.copy()
			self.sat3 = sat3.copy()
		else:
			lisa = LISAutils.LISA(initialConstellationPhase, initialLongitude)
			self.sat1 = State(lisa.sat(0).pos(), lisa.sat(0).vel())
			self.sat2 = State(lisa.sat(1).pos(), lisa.sat(1).vel())
			self.sat3 = State(lisa.sat(2).pos(), lisa.sat(2).vel())

	def copy(self):
		return LISA_State(self.sat1, self.sat2, self.sat3)

	def __add__(self, other):
		return LISA_State(self.sat1 + other.sat1,
			self.sat2 + other.sat2,
			self.sat3 + other.sat3)

	def __mul__(self, other):
		return LISA_State(self.sat1 * other,
			self.sat2 * other,
			self.sat3 * other)

	def __rmul__(self, other):
		return self * other

	def __truediv__(self, other):
		return self * (1 / other)

	def __neg__(self):
		return self * -1

	def __sub__(self, other):
		return self + (-other)


def solarForce(pos):
	posSqr = sum(pos * pos)
	return -pos * mu / posSqr**(1.5)

def bodyAccel(pos, time):
	otherPos = tools.linInterp(ceresPositions, 0, dt, time)
	displacement = pos - otherPos
	posSqr = sum(displacement**2)
	specMu = ceresMass * G
	return -displacement * specMu / posSqr**(1.5)

def solarDifference(pos, deltaPos):
	#returns solarForce(pos + deltaPos) - solarForce(pos)
	#up to second order in epsilon (defined below)
	posSqr = sum(pos**2)
	epsilon = 2 * sum(pos * deltaPos) / posSqr + sum(deltaPos**2) / posSqr
	# alpha = 3 / 2 * epsilon
	alpha = 3 / 2 * epsilon - 15 / 8 * epsilon**2
	coef = mu / posSqr**1.5
	output = pos * coef * alpha
	output -= deltaPos * coef * (1 - alpha)
	return output

def stateDerivative(state, time):
	return State(state.vel, solarForce(state.pos),
		state.dVel, bodyAccel(state.pos + state.dPos, time)
		+ solarDifference(state.pos, state.dPos))

def lisaStateDerivative(lisaState, time):
	return LISA_State(stateDerivative(lisaState.sat1, time),
		stateDerivative(lisaState.sat2, time),
		stateDerivative(lisaState.sat3, time))


def magSqr(inputVector):
	return sum(inputVector * inputVector)

def stateError(state, deltaState, dt):
	e1 = (magSqr(deltaState.pos) + dt**2 * magSqr(deltaState.vel)) /\
		(magSqr(state.pos) + dt**2 * magSqr(state.vel))
	e2 = (magSqr(deltaState.dPos) + dt**2 * magSqr(deltaState.dVel)) /\
		(magSqr(state.dPos) + dt**2 * magSqr(state.dVel) + 1)
	return epsilon / np.sqrt(e1 + e2)

def lisaStateError(lisaState, deltaState, dt):
	return min(stateError(lisaState.sat1, deltaState.sat1, dt),
		stateError(lisaState.sat2, deltaState.sat2, dt),
		stateError(lisaState.sat3, deltaState.sat3, dt))



class basicState:
	def __init__(self, pos, vel):
		self.pos = pos
		self.vel = vel

	def copy(self):
		return basicState(self.pos.copy(), self.vel.copy())

	def __add__(self, other):
		return basicState(self.pos + other.pos,
			self.vel + other.vel)

	def __mul__(self, other):
		return basicState(self.pos * other,
			self.vel * other)

	def __rmul__(self, other):
		return self * other

	def __truediv__(self, other):
		return self * (1 / other)

	def __neg__(self):
		return self * -1

	def __sub__(self, other):
		return self + (-other)

def twoBodyStateDerivative(currentState, time):
	pos = currentState.pos
	posSqr = sum(pos * pos)
	return basicState(currentState.vel, -pos * mu / posSqr**(1.5))

def withCeresStateDerivative(currentState, time):
	otherPos = tools.linInterp(ceresPositions, 0, dt, time)
	displacement = currentState.pos - otherPos
	posSqr = sum(displacement**2)
	specMu = ceresMass * G
	ceresAccel = -displacement * specMu / posSqr**(1.5)
	pos = currentState.pos
	posSqr = sum(pos * pos)
	solarAccel = -pos * mu / posSqr**(1.5)
	return basicState(currentState.vel, ceresAccel + solarAccel)

def basicStateError(currentState, deltaState, dt):
	e1 = (magSqr(deltaState.pos) + dt**2 * magSqr(deltaState.vel)) /\
		(magSqr(currentState.pos) + dt**2 * magSqr(currentState.vel))
	return epsilon / np.sqrt(e1)

class LISA_BasicState:
	def __init__(self, initialConstellationPhase, initialLongitude, sat3 = None):
		if sat3 != None:
			self.sat1 = initialConstellationPhase.copy()
			self.sat2 = initialLongitude.copy()
			self.sat3 = sat3.copy()
		else:
			lisa = LISAutils.LISA(initialConstellationPhase, initialLongitude)
			self.sat1 = basicState(lisa.sat(0).pos(), lisa.sat(0).vel())
			self.sat2 = basicState(lisa.sat(1).pos(), lisa.sat(1).vel())
			self.sat3 = basicState(lisa.sat(2).pos(), lisa.sat(2).vel())

	def copy(self):
		return LISA_BasicState(self.sat1.copy(), self.sat2.copy(), self.sat3.copy())

	def __add__(self, other):
		return LISA_BasicState(self.sat1 + other.sat1,
			self.sat2 + other.sat2,
			self.sat3 + other.sat3)

	def __mul__(self, other):
		return LISA_BasicState(self.sat1 * other,
			self.sat2 * other,
			self.sat3 * other)

	def __rmul__(self, other):
		return self * other

	def __truediv__(self, other):
		return self * (1 / other)

	def __neg__(self):
		return self * -1

	def __sub__(self, other):
		return self + (-other)

def LISA_twoBodyStateDerivative(currentState, time):
	return LISA_BasicState(twoBodyStateDerivative(currentState.sat1, time),
		twoBodyStateDerivative(currentState.sat2, time),
		twoBodyStateDerivative(currentState.sat3, time))

def LISA_withCeresStateDerivative(currentState, time):
	return LISA_BasicState(withCeresStateDerivative(currentState.sat1, time),
		withCeresStateDerivative(currentState.sat2, time),
		withCeresStateDerivative(currentState.sat3, time))

def LISA_basicStateError(currentState, deltaState, dt):
	return min(basicStateError(currentState.sat1, deltaState.sat1, dt),
		basicStateError(currentState.sat2, deltaState.sat2, dt),
		basicStateError(currentState.sat3, deltaState.sat3, dt))


class DualState:
	def __init__(self, initialConstellationPhase, initialLongitude, sat3 = None):
		if type(initialLongitude) == LISA_BasicState:
			self.s1 = initialConstellationPhase
			self.s2 = initialLongitude
		else:
			self.s1 = LISA_BasicState(initialConstellationPhase, initialLongitude)
			self.s2 = self.s1.copy()

	def copy(self):
		return DualState(self.s1.copy(), self.s2.copy())

	def __add__(self, other):
		return DualState(self.s1 + other.s1,
			self.s2 + other.s2)

	def __mul__(self, other):
		return DualState(self.s1 * other,
			self.s2 * other)

	def __rmul__(self, other):
		return self * other

	def __truediv__(self, other):
		return self * (1 / other)

	def __neg__(self):
		return self * -1

	def __sub__(self, other):
		return self + (-other)

def dualStateDerivative(currentState, time):
	return DualState(LISA_twoBodyStateDerivative(currentState.s1, time),
		LISA_withCeresStateDerivative(currentState.s2, time))

def dualStateError(currentState, deltaState, dt):
	return min(LISA_basicStateError(currentState.s1, deltaState.s1, dt),
		LISA_basicStateError(currentState.s2, deltaState.s2, dt))






def performSimulation(initialConstellationPhase, startYear, totalTime, inputDT):
	global ceresPositions, dt
	times = []
	ceresPositions = []
	dt = inputDT

	startEpoch = LISAutils.yearToEpochMJD(startYear)
	# print(startEpoch)

	earth = LISAutils.orbit(math.radians(200.7), 0.01671, math.radians(0),
		math.radians(-11.261), math.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58324.75, name = 'earth')
	ceres = LISAutils.orbit(math.radians(352.23), 0.0755, math.radians(10.5935),
		math.radians(80.3099), math.radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = ceresMass)

	times = np.arange(0, totalTime, dt)
	for t in times:
		ceres.setTime(t)
		ceresPositions.append(ceres.pos())
	ceresPositions = np.array(ceresPositions)
	ceres.setTime(0)

	# times, tempStates = rk4(DualState(initialConstellationPhase, np.pi * 1/6 + np.pi * 2 / 3),
	# 	dualStateDerivative, dualStateError, 0, totalTime, dt)
	times, tempStates = rk4(DualState(initialConstellationPhase, np.arctan2(earth.y(), earth.x())),
		dualStateDerivative, dualStateError, 0, totalTime, dt)

	stacks = [[],[],[]]
	for state in tempStates:
		p1 = state.s1.sat1.pos
		p2 = state.s1.sat2.pos
		p3 = state.s1.sat3.pos
		p1p = state.s2.sat1.pos
		p2p = state.s2.sat2.pos
		p3p = state.s2.sat3.pos
		v1 = state.s1.sat1.vel
		v2 = state.s1.sat2.vel
		v3 = state.s1.sat3.vel
		v1p = state.s2.sat1.vel
		v2p = state.s2.sat2.vel
		v3p = state.s2.sat3.vel
		stacks[0].append(State(p1, v1, p1p-p1, v1p-v1))
		stacks[1].append(State(p2, v2, p2p-p2, v2p-v2))
		stacks[2].append(State(p3, v3, p3p-p3, v3p-v3))

	newTimes = []
	for t in times:
		newTimes.append(t / year + startYear)
	newTimes = np.array(newTimes)

	return stacks, times, ceresPositions, newTimes, dt


def perturbativePerformSimulation(initialConstellationPhase, startYear, totalTime, inputDT):
	global ceresPositions, dt
	times = []
	ceresPositions = []
	dt = inputDT

	startEpoch = LISAutils.yearToEpochMJD(startYear)

	earth = LISAutils.orbit(math.radians(200.7), 0.01671, math.radians(0),
		math.radians(-11.261), math.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58324.75, name = 'earth')
	ceres = LISAutils.orbit(math.radians(352.23), 0.0755, math.radians(10.5935),
		math.radians(80.3099), math.radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = ceresMass)

	times = np.arange(0, totalTime, dt)
	for t in times:
		ceres.setTime(t)
		ceresPositions.append(ceres.pos())
	ceresPositions = np.array(ceresPositions)
	ceres.setTime(0)

	# times, tempStates = rk4(LISA_State(initialConstellationPhase, np.pi * 1/6 + np.pi * 2 / 3),
	# 	lisaStateDerivative, lisaStateError, 0, totalTime, dt)
	temp = LISA_State(initialConstellationPhase, np.arctan2(earth.y(), earth.x()))
	sats = []
	sats.append(temp.sat1)
	sats.append(temp.sat2)
	sats.append(temp.sat3)
	for sat in sats:
		pos = sat.pos
		cPos = ceresPositions[0]
		dPos = pos - cPos
		dist = np.sqrt(sum(dPos * dPos))
		print(ceres.mass * 6.67e-11 / dist)

	times, tempStates = rk4(LISA_State(initialConstellationPhase, np.arctan2(earth.y(), earth.x())),
		lisaStateDerivative, lisaStateError, 0, totalTime, dt)

	stacks = [[],[],[]]
	for state in tempStates:
		stacks[0].append(state.sat1)
		stacks[1].append(state.sat2)
		stacks[2].append(state.sat3)

	newTimes = []
	for t in times:
		newTimes.append(t / year + startYear)
	newTimes = np.array(newTimes)

	return stacks, times, ceresPositions, newTimes, dt




	# stacks = []
	# for i in range(3):
	# 	initialState = State(lisa.sat(i).pos(), lisa.sat(i).vel())
	# 	times, stack = rk4(initialState, stateDerivative, stateError, 0, totalTime, dt)