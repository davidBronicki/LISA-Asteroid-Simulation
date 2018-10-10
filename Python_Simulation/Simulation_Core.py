import numpy as np
import LISAutils
import math

G = 6.67408e-11
mu = G * 1.989e30
year = 365 * 24 * 3600
ceresMass = 8.958e20
ceresPositions = []
dt = 0

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


def performSimulation(initialConstellationPhase, startYear, totalTime, inputDT):
	global ceresPositions, dt
	times = []
	ceresPositions = []
	newTimes = []
	dt = inputDT

	startEpoch = LISAutils.yearToEpochMJD(startYear)

	earth = LISAutils.orbit(math.radians(200.7), 0.01671, math.radians(0),
		math.radians(-11.261), math.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58324.75, name = 'earth')
	lisa = LISAutils.LISA(initialConstellationPhase,
		np.arctan2(earth.y(), earth.x()))
		# np.pi * 1/6 + np.pi * 2 / 3)
	ceres = LISAutils.orbit(math.radians(352.23), 0.0755, math.radians(10.5935),
		math.radians(80.3099), math.radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
		simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = ceresMass)

	times = np.arange(0, totalTime, dt)
	for t in times:
		ceres.setTime(t)
		ceresPositions.append(ceres.pos())
		newTimes.append(t / year + startYear)
	newTimes = np.array(newTimes)
	ceresPositions = np.array(ceresPositions)
	ceres.setTime(0)

	stacks = []
	for i in range(3):
		initialState = State(lisa.sat(i).pos(), lisa.sat(i).vel())
		stacks.append(rk4(initialState, stateDerivative, 0, totalTime, dt))

	return stacks, times, ceresPositions, newTimes, dt