import numpy as np
import LISAutils
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# import matplotlib.colors as cl
import my_package.util as tools
from math import *

from Simulation_Core import *
from Data_Extraction_Methods import *


stacks, times, ceresPositions, newTimes, dt = performSimulation(np.pi / 6 * 0, 2029, 6*year, 1e5)

# relAngle = relativeAngleToCeres(stacks)
# normAngleToCeres = normAngleFromArmsToCeres(stacks)
# diff = get_rHatPhiHatGuidingCenterD_Pos(stacks)
# diff = getGuidingCenterD_Pos(stacks)
armDif = getD_ArmLengths(stacks)
deltaDistToGuide = getD_DistToGuide(stacks)
# armAccels = getArmAccels(stacks, times)
twoBodyPerturbAccels = getPerturbativeAccelValues(stacks)
# armLengths = getArmLengths(stacks)
# angleDif = getD_ArmAngles(stacks)
# armAccelDotVel = getArmAccelDotVel(stacks)



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

plt.figure(1)
plt.plot(newTimes, twoBodyPerturbAccels[0])
plt.plot(newTimes, twoBodyPerturbAccels[1])
plt.plot(newTimes, twoBodyPerturbAccels[2])
plt.grid()
plt.title('perturbative acceleration')
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
# plt.figure(4)
# plt.plot(newTimes, armDif[0])
# plt.plot(newTimes, armDif[1])
# plt.plot(newTimes, armDif[2])
# plt.grid()
# plt.title('arm perturbation')
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


plt.figure(1)
plt.plot(newTimes, deltaDistToGuide[0], label = "Satellite 1")
plt.plot(newTimes, deltaDistToGuide[1], label = "Satellite 2")
plt.plot(newTimes, deltaDistToGuide[2], label = "Satellite 3")
plt.grid()
plt.title('Perturbation of the Distance of LISA Satellites\nto the Guiding Center (Constellation Angle 120)')
plt.xlabel('time (year)')
plt.ylabel('displacement (meters)')
# plt.ylim([-18, 18])
plt.legend()


# plt.figure(1)
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
plt.show()