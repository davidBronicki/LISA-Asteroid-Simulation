import numpy as np
import matplotlib.pyplot as plt
import my_module.util as tools
import dataHandler as dh

handler = dh.dataHandler('such and such file')

times = handler.getTimes()
t0 = times[0]
dt = times[1] - times[0]
acc1 = tools.interpFunct(handler.getEclipticAcc(1), t0, dt)
acc2 = tools.interpFunct(handler.getEclipticAcc(2), t0, dt)
acc3 = tools.interpFunct(handler.getEclipticAcc(3), t0, dt)
accMat = lambda t: np.array([acc1(t), acc2(t), acc3(t)])

au = 1.5e11
G = 6.674e-11
solarM = 1.989e30
mu = G * solarM

def centeralMatForce(posMat):
	output = []
	for pos in posMat:
		output.append(-r * mu / np.sqrt(sum(pos * pos))**3)
	return np.array(output)

class LISA:
	def __init__(self, initialEarthAngle, initialOrientationAngle):
		pass

	def update(self, t, dt):
		self.perturbVelMat += dt * (centeralMatForce(
			self.perturbPosMat) + accMat(t))
		self.perturbPosMat += dt * self.perturbVelMat
		self.perturbMatList.append(self.perturbPosMat.copy())

		self.plainVelMat += dt * centeralMatForce(
			self.plainPosMat)
		self.plainPosMat += dt * self.plainVelMat
		self.plainMatList.append(self.plainPosMat.copy())



#####################################################
#######build initial LISA setup in state space#######
#####################################################


lisa = LISA(-20, 0)

dt = 1000
timeList = np.arange(0, 1000 * 24 * 3600, dt)

for t in timeList:
	lisa.update(t, dt)
