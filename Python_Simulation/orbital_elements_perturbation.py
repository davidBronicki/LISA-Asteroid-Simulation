import numpy as np
import LISAutils
# import my_package.util as tools
# import time as timePackage
from Data_Extraction_Methods import *
from matplotlib import pyplot as plt

G = 6.67408e-11
mu = G * 1.989e30
year = 365 * 24 * 3600
ceresMass = 8.958e20
epsilon = 1e-10


degreesOfPhase = 0

numberOfYear = 5
dt = 1e5
initYear = 2029
length = year*numberOfYear
constellationPhase = -np.pi / 3 + np.radians(degreesOfPhase)
# initPhase = np.pi / 6

# stacks, times, ceresPositions, newTimes, dt = performSimulation(initPhase, initYear, length, dt)


startEpoch = LISAutils.yearToEpochMJD(initYear)
J2000 = 51544

earth = LISAutils.orbit(np.radians(200.7), 0.01671, np.radians(0),
	np.radians(-11.261), np.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58324.75, name = 'Earth', mass = 5.972e24)
ceres = LISAutils.orbit(np.radians(352.23), 0.0755, np.radians(10.5935),
	np.radians(80.3099), np.radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58200, name = 'Ceres', mass = ceresMass)

mars = LISAutils.orbit(np.radians(355.45), 0.0934, np.radians(1.85),
	np.radians(49.58), np.radians(336.04), 1.524 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = J2000, name = 'Mars', mass = 6.39e23)
venus = LISAutils.orbit(np.radians(181.97973), 0.00677323, np.radians(3.39471),
	np.radians(76.68069), np.radians(131.53298), 0.7233 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = J2000, name = 'Venus', mass = 4.867e24)

jupiter = LISAutils.orbit(np.radians(34.4), 0.048, np.radians(1.305),
	np.radians(100.556), np.radians(14.754), 5.2 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = J2000, name = 'Jupiter', mass = 1.898e27)
saturn = LISAutils.orbit(np.radians(49.94432), 0.05415060, np.radians(2.48446),
	np.radians(113.71504), np.radians(92.43194), 9.54 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = J2000, name = 'Saturn', mass = 5.683e26)


def performSimulation(perturbingBody, flipLISA = False):
	lisa_unperturbed = LISAutils.LISA(constellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9, flipLISA)
	lisa_perturbed = LISAutils.LISA(constellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9, flipLISA)

	stacks = [[],[],[]]
	times = []
	perturbingBodyPositions = []
	# elements = [[],[],[],[],[],[]]#anom, semi-major, ecc, ang to par, inc, RAAN

	for t in np.arange(0, length, dt):
		lisa_perturbed.iterateTime(dt)
		lisa_unperturbed.iterateTime(dt)
		# lisa_perturbed.setTime(t)
		# lisa_unperturbed.setTime(t)
		perturbingBody.iterateTime(dt)

		times.append(t)
		perturbingBodyPositions.append(perturbingBody.pos())

			# elements[0].append(lisa_perturbed.sats[0].f-lisa_unperturbed.sats[0].f)
			# elements[1].append(
			# 	(lisa_perturbed.sats[0].a-lisa_unperturbed.sats[0].a)
			# 	/lisa_unperturbed.sats[0].a)
			# elements[2].append(
			# 	(lisa_perturbed.sats[0].e-lisa_unperturbed.sats[0].e)
			# 	/lisa_unperturbed.sats[0].e)
			# elements[3].append(
			# 	(lisa_perturbed.sats[0].omega-lisa_unperturbed.sats[0].omega)
			# 	/lisa_unperturbed.sats[0].omega)
			# elements[4].append(
			# 	(lisa_perturbed.sats[0].i-lisa_unperturbed.sats[0].i)
			# 	/lisa_unperturbed.sats[0].i)
			# elements[5].append(
			# 	(lisa_perturbed.sats[0].RAAN-lisa_unperturbed.sats[0].RAAN)
			# 	/lisa_unperturbed.sats[0].RAAN)

		dp = lisa_unperturbed.pos()
		dpp = lisa_perturbed.pos()
		dv = lisa_unperturbed.vel()
		dvp = lisa_perturbed.vel()

		stacks[0].append(State(dp[0],dv[0],dpp[0]-dp[0],dvp[0]-dv[0]))
		stacks[1].append(State(dp[1],dv[1],dpp[1]-dp[1],dvp[1]-dv[1]))
		stacks[2].append(State(dp[2],dv[2],dpp[2]-dp[2],dvp[2]-dv[2]))

		# if t < 2*np.pi*1e7:
		lisa_perturbed.beginPerturb()
		lisa_perturbed.perturbFrom(perturbingBody, dt)
		lisa_perturbed.endPerturb()

	newTimes = []
	for t in times:
		newTimes.append(t / year + initYear)
	newTimes = np.array(newTimes)

	return (stacks, newTimes, perturbingBodyPositions)

def performStaticSimulation(distanceToSun):
	lisa_unperturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)
	lisa_perturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)

	perturbingBody = LISAutils.orbit(0, 0, np.radians(0.1),
		0, 0, distanceToSun * LISAutils.astroUnit, anomType = 'trueAnom',
		simTime = startEpoch, paramTime = startEpoch, name = 'Static Mass', mass = 5e26)

	stacks = [[],[],[]]
	times = []

	for t in np.arange(0, length, dt):
		lisa_perturbed.iterateTime(dt)
		lisa_unperturbed.iterateTime(dt)
		perturbingBody.iterateTime(dt)

		times.append(t)

		dp = lisa_unperturbed.pos()
		dpp = lisa_perturbed.pos()
		dv = lisa_unperturbed.vel()
		dvp = lisa_perturbed.vel()

		stacks[0].append(State(dp[0],dv[0],dpp[0]-dp[0],dvp[0]-dv[0]))
		stacks[1].append(State(dp[1],dv[1],dpp[1]-dp[1],dvp[1]-dv[1]))
		stacks[2].append(State(dp[2],dv[2],dpp[2]-dp[2],dvp[2]-dv[2]))

		lisa_perturbed.beginPerturb()
		lisa_perturbed.perturbFrom(perturbingBody, dt)
		lisa_perturbed.endPerturb()

	newTimes = []
	for t in times:
		newTimes.append(t / year + initYear)
	newTimes = np.array(newTimes)

	return (stacks, newTimes)

def performImpulseSimulation(distanceToSun):
	lisa_unperturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)
	lisa_perturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)

	perturbingBody = LISAutils.orbit(0, 0, np.radians(0.1),
		0, 0, distanceToSun * LISAutils.astroUnit, anomType = 'trueAnom',
		simTime = startEpoch, paramTime = startEpoch, name = 'Static Mass', mass = 5e27)

	stacks = [[],[],[]]
	times = []

	for t in np.arange(0, length, dt):
		lisa_perturbed.iterateTime(dt)
		lisa_unperturbed.iterateTime(dt)
		perturbingBody.iterateTime(dt)

		times.append(t)

		dp = lisa_unperturbed.pos()
		dpp = lisa_perturbed.pos()
		dv = lisa_unperturbed.vel()
		dvp = lisa_perturbed.vel()

		stacks[0].append(State(dp[0],dv[0],dpp[0]-dp[0],dvp[0]-dv[0]))
		stacks[1].append(State(dp[1],dv[1],dpp[1]-dp[1],dvp[1]-dv[1]))
		stacks[2].append(State(dp[2],dv[2],dpp[2]-dp[2],dvp[2]-dv[2]))

		if t == 0:
			lisa_perturbed.beginPerturb()
			lisa_perturbed.perturbFrom(perturbingBody, dt)
			lisa_perturbed.endPerturb()

	newTimes = []
	for t in times:
		newTimes.append(t / year + initYear)
	newTimes = np.array(newTimes)

	return (stacks, newTimes)

def performManyBodySimulations(perturbingBodies):
	lisa_unperturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)
	lisa_perturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()) - np.pi / 9)

	stacks = [[],[],[]]
	times = []
	# elements = [[],[],[],[],[],[]]#anom, semi-major, ecc, ang to par, inc, RAAN

	for t in np.arange(0, length, dt):
		lisa_perturbed.iterateTime(dt)
		lisa_unperturbed.iterateTime(dt)

		times.append(t)

		dp = lisa_unperturbed.pos()
		dpp = lisa_perturbed.pos()
		dv = lisa_unperturbed.vel()
		dvp = lisa_perturbed.vel()

		stacks[0].append(State(dp[0],dv[0],dpp[0]-dp[0],dvp[0]-dv[0]))
		stacks[1].append(State(dp[1],dv[1],dpp[1]-dp[1],dvp[1]-dv[1]))
		stacks[2].append(State(dp[2],dv[2],dpp[2]-dp[2],dvp[2]-dv[2]))

		for body in perturbingBodies:
			body.iterateTime(dt)
			lisa_perturbed.beginPerturb()
			lisa_perturbed.perturbFrom(body, dt)
			lisa_perturbed.endPerturb()

	newTimes = []
	for t in times:
		newTimes.append(t / year + initYear)
	newTimes = np.array(newTimes)

	return (stacks, newTimes)

def staticBodyRunner():
	stacks, newTimes = performStaticSimulation(3)

	armLengths = getArmLengths(stacks)
	armDif = getD_ArmLengths(stacks)

	plt.figure(0)
	# plt.plot(newTimes, armLengths[0]/1e9, label = 'arm 1')
	plt.plot(newTimes, armLengths[1]/1e9, label = 'arm 2')
	# plt.plot(newTimes, armLengths[2]/1e9, label = 'arm 3')
	plt.grid()
	plt.legend()
	plt.xlabel('time (years)')
	plt.ylabel('length (Gm)')
	plt.title('Length of Arm 2 in Presence of Static Body')

	plt.figure(2)
	plt.plot(newTimes, armDif[0]/1e6, label = 'arm 1 displacement')
	plt.plot(newTimes, armDif[1]/1e6, label = 'arm 2 displacement')
	plt.plot(newTimes, armDif[2]/1e6, label = 'arm 3 displacement')
	plt.ylabel('displacement (Mm)')
	plt.grid()
	plt.xlabel('time (years)')
	plt.title('Perturbation of Arm Lengths due to Static Body\n(5e27 kg mass at 3 au)')
	plt.legend()
	plt.ylim([-5,5])

	plt.show()

def impulseRunner():
	stacks, newTimes = performImpulseSimulation(3)

	armLengths = getArmLengths(stacks)
	armDif = getD_ArmLengths(stacks)

	plt.figure(0)
	# plt.plot(newTimes, armLengths[0]/1e9, label = 'arm 1')
	plt.plot(newTimes, armLengths[1]/1e9, label = 'arm 2')
	# plt.plot(newTimes, armLengths[2]/1e9, label = 'arm 3')
	plt.grid()
	plt.legend()
	plt.xlabel('time (years)')
	plt.ylabel('length (Gm)')
	plt.title('Length of Arm 2 in Presence of Static Body')

	plt.figure(2)
	plt.plot(newTimes, armDif[0]/1e6, label = 'arm 1 displacement')
	plt.plot(newTimes, armDif[1]/1e6, label = 'arm 2 displacement')
	plt.plot(newTimes, armDif[2]/1e6, label = 'arm 3 displacement')
	plt.ylabel('displacement (Mm)')
	plt.grid()
	plt.xlabel('time (years)')
	plt.title('Perturbation of Arm Lengths due to Momentary Impulse\n(5e27 kg mass at 3 au applied for 1e5 s)')
	plt.legend()
	# plt.ylim([-5,5])

	plt.show()

def planetRunner():
	earthStack, newTimes, throwAway = performSimulation(earth)
	venusStack, newTimes, throwAway = performSimulation(venus)
	jupiterStack, newTimes, throwAway = performSimulation(jupiter)
	saturnStack, newTimes, throwAway = performSimulation(saturn)
	marsStack, newTimes, throwAway = performSimulation(mars)
	# earthDif = getD_ArmLengths(earthStack)
	# venusDif = getD_ArmLengths(venusStack)
	# jupiterDif = getD_ArmLengths(jupiterStack)
	# saturnDif = getD_ArmLengths(saturnStack)
	# marsDif = getD_ArmLengths(marsStack)

	# area = constellationArea(stacks)
	dAreaEarth = getPerburbationOfArea(earthStack)
	dAreaVenus = getPerburbationOfArea(venusStack)
	dAreaJupiter = getPerburbationOfArea(jupiterStack)
	dAreaSaturn = getPerburbationOfArea(saturnStack)
	dAreaMars = getPerburbationOfArea(marsStack)

	# plt.figure(0)
	# plt.plot(newTimes, earthDif[1]/1e6, label = 'from Earth')
	# plt.plot(newTimes, venusDif[1]/1e6, label = 'from Venus')
	# plt.plot(newTimes, jupiterDif[1]/1e6, label = 'from Jupiter')
	# plt.plot(newTimes, saturnDif[1]/1e6, label = 'from Saturn')
	# plt.plot(newTimes, marsDif[1]/1e6, label = 'from Mars')
	# plt.grid()
	# plt.legend()
	# plt.xlabel('time (years)')
	# plt.ylabel('displacement (Mm)')
	# plt.title('Perturbation of Arm 2 from Planets')

	plt.figure(1)
	plt.plot(newTimes, dAreaEarth/1e15, label = 'from Earth')
	plt.plot(newTimes, dAreaVenus/1e15, label = 'from Venus')
	plt.plot(newTimes, dAreaJupiter/1e15, label = 'from Jupiter')
	plt.plot(newTimes, dAreaSaturn/1e15, label = 'from Saturn')
	plt.plot(newTimes, dAreaMars/1e15, label = 'from Mars')
	plt.grid()
	plt.legend()
	plt.ylabel('$\delta$area ($10^{15}\, m^2$)')
	plt.xlabel('time (years)')
	plt.title('Perturbation of the Enclosed Area of LISA\nDue to Various Planets')

	plt.show()
'''
do with different configurations
screen capture and send to shane
'''

# planetRunner()

def planetFullEffect():
	stacks, newTimes = performManyBodySimulations([earth, venus, jupiter, saturn, mars])

	armLengths = getArmLengths(stacks)
	armDif = getD_ArmLengths(stacks)

	plt.figure(0)
	# plt.plot(newTimes, armLengths[0]/1e9, label = 'arm 1')
	plt.plot(newTimes, armLengths[1]/1e9, label = 'arm 2')
	# plt.plot(newTimes, armLengths[2]/1e9, label = 'arm 3')
	plt.grid()
	plt.legend()
	plt.xlabel('time (years)')
	plt.ylabel('length (Gm)')
	plt.title('Length of Arm 2 in Presence of Planets')

	plt.figure(2)
	# plt.plot(newTimes, armDif[0]/1e6, label = 'arm 1 displacement')
	plt.plot(newTimes, armDif[1]/1e6, label = 'arm 2 displacement')
	# plt.plot(newTimes, armDif[2]/1e6, label = 'arm 3 displacement')
	plt.ylabel('displacement (Mm)')
	plt.grid()
	plt.xlabel('time (years)')
	plt.title('Perturbation of Arm 2 due to Planets')
	# plt.legend()
	plt.ylim([-22,22])

	plt.show()

def singleRunner():
	stacks, newTimes, perturbingBodyPositions = performSimulation(ceres)


	# relAngle = relativeAngleToCeres(stacks)
	# normAngleToCeres = normAngleFromArmsToCeres(stacks, times, dt, ceresPositions)
	# distanceToPerturbingBody = distToCeres(stacks, perturbingBodyPositions)
	# diff = get_rHatPhiHatGuidingCenterD_Pos(stacks)
	# diff = getGuidingCenterD_Pos(stacks)
	# armDif = getD_ArmLengths(stacks)
	# deltaDistToGuide = getD_DistToGuide(stacks)
	# armAccels = getArmAccels(stacks, times, ceresPositions)
	# twoBodyPerturbAccels = getPerturbativeAccelValues(stacks)
	# armLengths = getArmLengths(stacks)
	# area = constellationArea(stacks)
	dArea = getPerburbationOfArea(stacks)

	# dArmLengthSum = sum(armDif)*2.5e9/(2*np.sqrt(3))


	# plt.figure(0)
	# plt.plot(newTimes, armLengths[0]/1e9, label = 'arm 1')
	# plt.plot(newTimes, armLengths[1]/1e9, label = 'arm 2')
	# plt.plot(newTimes, armLengths[2]/1e9, label = 'arm 3')
	# plt.grid()
	# plt.legend()
	# plt.xlabel('time (years)')
	# plt.ylabel('length (Gm)')
	# plt.title('Arm Lengths')

	# plt.figure(1)
	# plt.plot(newTimes, area/1e18)
	# plt.xticks([2029,2030,2031])
	# plt.grid()
	# plt.ylabel("area ($10^{18}\, m^2$)")
	# plt.xlabel('time (years)')
	# plt.title("Enclosed Area of LISA Constellation\nWith Perturbations from Ceres")

	# plt.figure(3)
	# plt.plot(newTimes, dArmLengthSum/1e9)
	# plt.grid()
	# plt.ylabel('$\delta$ area ($10^9\, m^2$)')
	# plt.xlabel('time (years)')
	# plt.title('Perturbation of the Enclosed Area of LISA\nDue to Ceres at '\
	# 	+ str(degreesOfPhase) + '$\degree$ (Via Sum of Arm Length Change)')
	# plt.ylim([-20, 20])

	plt.figure(4)
	plt.plot(newTimes, dArea/1e9)
	plt.grid()
	plt.ylabel('$\delta$ area ($10^9\, m^2$)')
	plt.xlabel('time (years)')
	plt.title('Perturbation of the Enclosed Area of LISA\nDue to Ceres at '\
		+ str(degreesOfPhase) + '$\degree$ (Via Area Graph Difference)')
	plt.ylim([-20, 20])


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
	# plt.plot(newTimes, armDif[0], label = 'arm 1 variation')
	# plt.plot(newTimes, armDif[1], label = 'arm 2 variation')
	# plt.plot(newTimes, armDif[2], label = 'arm 3 variation')
	# plt.ylabel('$\delta$ Length (m)')
	# # plt.plot(newTimes, armDif[0]/1e3, label = 'arm 1 displacement')
	# # plt.plot(newTimes, armDif[1]/1e3, label = 'arm 2 displacement')
	# # plt.plot(newTimes, armDif[2]/1e3, label = 'arm 3 displacement')
	# # plt.ylabel('displacement (km)')
	# # plt.plot(newTimes, armDif[0]/1e6, label = 'arm 1 displacement')
	# # plt.plot(newTimes, armDif[1]/1e6, label = 'arm 2 displacement')
	# # plt.plot(newTimes, armDif[2]/1e6, label = 'arm 3 displacement')
	# # plt.ylabel('displacement (Mm)')
	# plt.grid()
	# plt.xlabel('time (years)')
	# # plt.title('Perturbation of LISA Arm Lengths due to Mars')#max 500 km
	# # plt.title('Perturbation of LISA Arm Lengths due to Jupiter')#max 2300 km
	# # plt.title('Perturbation of LISA Arm Lengths due to Earth')#max 23000 km
	# # plt.title('Perturbation of LISA Arm Lengths due to Venus')#max 9000 km
	# # plt.title('Perturbation of LISA Arm Lengths due to Saturn')#maxx 350 km
	# plt.title('Perturbation of LISA Arm Lengths due to Ceres\nat '\
	# 	+ str(degreesOfPhase) + '$\degree$ Constellation Phase Angle')
	# plt.legend()
	# plt.ylim([-20, 20])

	# plt.figure(3)
	# plt.plot(newTimes, distanceToPerturbingBody/LISAutils.astroUnit)
	# plt.grid()
	# plt.xlabel("time (years)")
	# plt.ylabel("distance (au)")
	# plt.title("distance from perturbing body")

	# plt.figure(3)
	# plt.plot(newTimes, elements[0])
	# plt.title('Perturbation of the True Anomaly')
	# plt.grid()

	# plt.figure(4)
	# plt.plot(newTimes, elements[1])
	# plt.title('Perturbation of the Semi-major Axis')
	# plt.grid()

	# plt.figure(5)
	# plt.plot(newTimes, elements[2], label = 'eccentricity')
	# plt.plot(newTimes, elements[3], label = 'argument of the pariapsis')
	# plt.plot(newTimes, elements[4], label = 'inclination')
	# plt.plot(newTimes, elements[5], label = 'right ascension of the ascending node')
	# plt.grid()
	# plt.title("Other Perturbations")
	# plt.legend()

	plt.show()

# singleRunner()
# planetRunner()

def areaFull():
	earthStack, newTimes, throwAway = performSimulation(earth)
	venusStack, newTimes, throwAway = performSimulation(venus)
	jupiterStack, newTimes, throwAway = performSimulation(jupiter)
	saturnStack, newTimes, throwAway = performSimulation(saturn)
	marsStack, newTimes, throwAway = performSimulation(mars)

	dAreaEarthRegular = getPerburbationOfArea(earthStack)
	dAreaVenusRegular = getPerburbationOfArea(venusStack)
	dAreaJupiterRegular = getPerburbationOfArea(jupiterStack)
	dAreaSaturnRegular = getPerburbationOfArea(saturnStack)
	dAreaMarsRegular = getPerburbationOfArea(marsStack)

	plt.figure(1)
	plt.plot(newTimes, dAreaEarthRegular/1e15, label = 'from Earth')
	plt.plot(newTimes, dAreaVenusRegular/1e15, label = 'from Venus')
	plt.plot(newTimes, dAreaJupiterRegular/1e15, label = 'from Jupiter')
	plt.plot(newTimes, dAreaSaturnRegular/1e15, label = 'from Saturn')
	plt.plot(newTimes, dAreaMarsRegular/1e15, label = 'from Mars')
	plt.grid()
	plt.legend()
	plt.ylabel('$\delta$area ($10^{15}\, m^2$)')
	plt.xlabel('time (years)')
	plt.title('Perturbation of the Enclosed Area of LISA\nDue to Various Planets with Ordinary Orientation')



	earthStack, newTimes, throwAway = performSimulation(earth, True)
	venusStack, newTimes, throwAway = performSimulation(venus, True)
	jupiterStack, newTimes, throwAway = performSimulation(jupiter, True)
	saturnStack, newTimes, throwAway = performSimulation(saturn, True)
	marsStack, newTimes, throwAway = performSimulation(mars, True)

	dAreaEarthFlipped = getPerburbationOfArea(earthStack)
	dAreaVenusFlipped = getPerburbationOfArea(venusStack)
	dAreaJupiterFlipped = getPerburbationOfArea(jupiterStack)
	dAreaSaturnFlipped = getPerburbationOfArea(saturnStack)
	dAreaMarsFlipped = getPerburbationOfArea(marsStack)

	plt.figure(2)
	plt.plot(newTimes, dAreaEarthFlipped/1e15, label = 'from Earth')
	plt.plot(newTimes, dAreaVenusFlipped/1e15, label = 'from Venus')
	plt.plot(newTimes, dAreaJupiterFlipped/1e15, label = 'from Jupiter')
	plt.plot(newTimes, dAreaSaturnRegular/1e15, label = 'from Saturn')
	plt.plot(newTimes, dAreaMarsRegular/1e15, label = 'from Mars')
	plt.grid()
	plt.legend()
	plt.ylabel('$\delta$area ($10^{15}\, m^2$)')
	plt.xlabel('time (years)')
	plt.title('Perturbation of the Enclosed Area of LISA\nDue to Various Planets with Flipped Orientation')



	ceresStacks, newTimes, throwAway = performSimulation(ceres)
	dAreaCeresRegular = getPerburbationOfArea(ceresStacks)
	ceresStacks, newTimes, throwAway = performSimulation(ceres, True)
	dAreaCeresFlipped = getPerburbationOfArea(ceresStacks)
	
	plt.figure(0)
	plt.plot(newTimes, dAreaCeresRegular/1e9, label = "Regular")
	plt.plot(newTimes, dAreaCeresFlipped/1e9, label = "Flipped")
	plt.grid()
	plt.legend()
	plt.ylabel('$\delta$ area ($10^9\, m^2$)')
	plt.xlabel('time (years)')
	plt.title('Perturbation of the Enclosed Area of LISA\nDue to Ceres with Different Parities')
	plt.ylim([-20, 20])

	plt.show()

areaFull()
