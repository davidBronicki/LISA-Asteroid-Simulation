import numpy as np
from matplotlib import pyplot as plt

def funct(x0, x, phase, linSlope, phaseRate):
	return linSlope * (x-x0) * np.sin(x * phaseRate + phase)

def display(x0, phase, phaseRate, m1, m2, m3):
	xList = np.arange(2029, 2035, 0.01)
	plt.plot(xList, funct(x0, xList, -phase, m2, phaseRate))
	plt.plot(xList, funct(x0, xList, -phase-np.pi/3, m3, phaseRate))
	plt.plot(xList, funct(x0, xList, -phase+np.pi/3, m1, phaseRate))
	plt.grid()
	plt.ylim([-20,20])
	plt.show()

display(2029.5, 2*np.pi*0.3, 2*np.pi, 4.2, 5.8, 1.8)
