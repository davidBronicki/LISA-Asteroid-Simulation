import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# import matplotlib.colors as cl
import my_package.util as tools
# from math import *
import math

from Simulation_Core import *
from Data_Extraction_Methods import *

simulationDataDirectory = "../../Build_Directory/"

fileName = "outputData.csv"

data = tools.floatParseCSV(simulationDataDirectory + fileName)

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

newTimes = times

# relAngle = relativeAngleToCeres(stacks)
diff = get_rHatPhiHatGuidingCenterD_Pos(stacks)
# diff = getGuidingCenterD_Pos(stacks)
armDif = getD_ArmLengths(stacks)
deltaDistToGuide = getD_DistToGuide(stacks)
twoBodyPerturbAccels = getPerturbativeAccelValues(stacks)
armLengths = getArmLengths(stacks)
# angleDif = getD_ArmAngles(stacks)
# armAccelDotVel = getArmAccelDotVel(stacks)


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


# plt.figure(2)
# # plt.plot(newTimes, relAngle, label = 'angle to Ceres', linestyle=':')
# plt.plot(newTimes, diff[0], label = 'r displacement')
# # plt.plot(newTimes, diff[1], label = 'phi displacement')
# # plt.plot(newTimes, diff[2], label = 'z displacement')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('displacement (meters)')
# plt.title('Perturbation of LISA Guiding Center due to Ceres\nin r-hat Direction')
# # plt.ylim([-180, 180])
# # plt.legend()

plt.figure(3)
plt.plot(newTimes, armDif[0], label = 'arm 1 displacement')
plt.plot(newTimes, armDif[1], label = 'arm 2 displacement')
plt.plot(newTimes, armDif[2], label = 'arm 3 displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title("Perturbation of LISA Arm Lengths due to 1k Asteroids")
# plt.title("Perturbation of LISA Arm Lengths due to Ceres")
# plt.title('Perturbation of LISA Arm Lengths due to Ceres\nat 0 Degree Initial Phase Angle')#\nwith 0 constellation phase angle')
plt.legend()
plt.ylim([-20, 20])

plt.figure(4)
plt.plot(newTimes, armLengths[0], label = 'arm 1 length')
plt.plot(newTimes, armLengths[1], label = 'arm 2 length')
plt.plot(newTimes, armLengths[2], label = 'arm 3 length')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title("LISA Arm Lengths")
plt.legend()

# plt.figure(3)
# plt.plot(newTimes, normAngleToCeres[0]*180/np.pi, label = 'arm 1 normal angle to Ceres')
# plt.plot(newTimes, normAngleToCeres[1]*180/np.pi, label = 'arm 2 normal angle to Ceres')
# plt.plot(newTimes, normAngleToCeres[2]*180/np.pi, label = 'arm 3 normal angle to Ceres')
# # plt.plot(newTimes, relAngle, color = 'black')# label = 'half angle from guiding center to Ceres', color = 'black')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('angle (degrees)')
# plt.title('Relative Angles Between Normal and Ceres')# to LISA Arm 1 and Ceres. Black line is angle from LISA Constellation to Ceres\nColor Mapped According to Constellation Angle on x-axis: red = 0 degrees, blue = 60 degrees')#\nwith 0 constellation phase angle')
# # plt.legend()
# plt.ylim([-90, 90])

# plt.figure(4)
# plt.plot(newTimes, 10**15*armAccels[0], label = 'arm 1 acceleration')
# plt.plot(newTimes, 10**15*armAccels[1], label = 'arm 2 acceleration')
# plt.plot(newTimes, 10**15*armAccels[2], label = 'arm 3 acceleration')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('acceleration (pN/kg)')
# plt.title('Perturbative Acceleration due to Ceres')
# plt.legend()
# plt.ylim([-12, 12])

# plt.figure(5)
# plt.plot(newTimes, angleDif[0] * 180 / np.pi * 3600, label = 'angle 1 displacement')
# plt.plot(newTimes, angleDif[1] * 180 / np.pi * 3600, label = 'angle 2 displacement')
# plt.plot(newTimes, angleDif[2] * 180 / np.pi * 3600, label = 'angle 3 displacement')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('angle (arcseconds)')
# plt.title('Perturbation of LISA Arm Angles due to Ceres')
# plt.legend()
# plt.ylim([-0.0015, 0.0015])

# plt.figure(6)
# plt.plot(newTimes, distanceToCeres / LISAutils.astroUnit)
# # plt.plot(newTimes, armAccels[1], label = 'arm 2 acceleration')
# # plt.plot(newTimes, armAccels[2], label = 'arm 3 acceleration')
# plt.grid()
# plt.xlabel('time (years)')
# plt.ylabel('distance (au)')
# plt.title('Distance from LISA to Ceres')
# plt.legend()
# # plt.ylim([-12, 12])


plt.show()

