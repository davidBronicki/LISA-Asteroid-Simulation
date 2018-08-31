import numpy as np
import LISAutils
import matplotlib.pyplot as plt

mu = 6.67408e-11 * 1.989e30

def eulerVelFirst(initPos, initVel, accelFunct, finalT, dt):
	timeList = np.arange(0, finalT, dt)
	pos = initPos.copy()
	posList = [initPos.copy()]
	vel = initVel.copy()
	for t in timeList:
		vel += accelFunct(pos) * dt
		pos += vel * dt
		posList.append(pos.copy())
	return posList

def eulerPosFirst(initPos, initVel, accelFunct, finalT, dt):
	timeList = np.arange(0, finalT, dt)
	pos = initPos.copy()
	posList = [initPos.copy()]
	vel = initVel.copy()
	for t in timeList:
		pos += vel * dt
		vel += accelFunct(pos) * dt
		posList.append(pos.copy())
	return posList

def rk4(initPos, initVel, accelFunct, finalT, dt):
	timeList = np.arange(0, finalT, dt)
	pos = initPos.copy()
	posList = [initPos.copy()]
	vel = initVel.copy()
	for t in timeList:
		dv1 = dt * accelFunct(pos)
		dr1 = dt * vel
		dv2 = dt * accelFunct(pos + dr1 / 2)
		dr2 = dt * (vel + dv1 / 2)
		dv3 = dt * accelFunct(pos + dr2 / 2)
		dr3 = dt * (vel + dv2 / 2)
		dv4 = dt * accelFunct(pos + dr3)
		dr4 = dt * (vel + dv3)
		vel += (dv1 + 2 * dv2 + 2 * dv3 + dv4) / 6
		pos += (dr1 + 2 * dr2 + 2 * dr3 + dr4) / 6
		posList.append(pos.copy())
	return posList

def orbitalAccelFunct(pos):
	distSqr = sum(pos * pos)
	return -pos * mu / distSqr**(1.5)

earthish = LISAutils.orbit(0, 0, 0, 0, 0, 1.5e11)

year = 365 * 24 * 3600
dt = 1e5

posListEulerVelFirst = eulerVelFirst(earthish.pos(), earthish.vel(), orbitalAccelFunct,
	year*5, dt)
posListEulerPosFirst = eulerPosFirst(earthish.pos(), earthish.vel(), orbitalAccelFunct,
	year*5, dt)
posListRK4 = rk4(earthish.pos(), earthish.vel(), orbitalAccelFunct,
	year*5, dt)
posListAccual = [earthish.pos()]
for t in np.arange(dt, year*5 + dt, dt):
	earthish.setTime(t)
	posListAccual.append(earthish.pos())

timeList = np.arange(0, year*5 + dt, dt)

velFirstDif = np.array(posListEulerVelFirst) - np.array(posListAccual)
posFirstDif = np.array(posListEulerPosFirst) - np.array(posListAccual)
rk4Dif = np.array(posListRK4) - np.array(posListAccual)

plt.plot(timeList, posFirstDif, label = 'pos first')
plt.plot(timeList, velFirstDif, label = 'vel first')
plt.plot(timeList, rk4Dif, label = 'RK4')
plt.legend()
plt.show()
