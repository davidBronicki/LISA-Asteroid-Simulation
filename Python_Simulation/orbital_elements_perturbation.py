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


numberOfYear = 5
dt = 1e5
initYear = 2029
length = year*numberOfYear
initialConstellationPhase = -np.pi / 3 + np.radians(80)
# initPhase = np.pi / 6

# stacks, times, ceresPositions, newTimes, dt = performSimulation(initPhase, initYear, length, dt)


startEpoch = LISAutils.yearToEpochMJD(initYear)

earth = LISAutils.orbit(np.radians(200.7), 0.01671, np.radians(0),
	np.radians(-11.261), np.radians(114.2078), LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58324.75, name = 'earth')
ceres = LISAutils.orbit(np.radians(352.23), 0.0755, np.radians(10.5935),
	np.radians(80.3099), np.radians(73.11534), 2.767 * LISAutils.astroUnit, anomType = 'meanAnom',
	simTime = startEpoch, paramTime = 58200, name = 'ceres', mass = ceresMass)


lisa_unperturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()))
lisa_perturbed = LISAutils.LISA(initialConstellationPhase, np.arctan2(earth.y(), earth.x()))



stacks = [[],[],[]]
times = []
ceresPositions = []
# elements = [[],[],[],[],[],[]]#anom, semi-major, ecc, ang to par, inc, RAAN

for t in np.arange(0, length, dt):
	lisa_perturbed.iterateTime(dt)
	lisa_unperturbed.iterateTime(dt)
	ceres.iterateTime(dt)

	times.append(t)
	ceresPositions.append(ceres.pos())

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

	lisa_perturbed.beginPerturb()
	lisa_perturbed.perturbFrom(ceres, dt)
	lisa_perturbed.endPerturb()

newTimes = []
for t in times:
	newTimes.append(t / year + initYear)
newTimes = np.array(newTimes)




posList = getPos(stacks)[0]
p_posList=posList + getD_Pos(stacks)[0]


accelSun = (posList[2]-2*posList[1]+posList[0])/(times[1]-times[0])**2
accel = (p_posList[2]-2*p_posList[1]+p_posList[0])/(times[1]-times[0])**2
print("Pos: " + str(posList[0]))
print("Pos relative to Ceres: " + str(posList[0]-ceresPositions[0]))
print("Acceleration from Ceres: " + str(accel-accelSun))





# relAngle = relativeAngleToCeres(stacks)
normAngleToCeres = normAngleFromArmsToCeres(stacks, times, dt, ceresPositions)
distanceToCeres = distToCeres(stacks, times, dt, ceresPositions)
diff = get_rHatPhiHatGuidingCenterD_Pos(stacks)
# diff = getGuidingCenterD_Pos(stacks)
armDif = getD_ArmLengths(stacks)
deltaDistToGuide = getD_DistToGuide(stacks)
# armAccels = getArmAccels(stacks, times, ceresPositions)
twoBodyPerturbAccels = getPerturbativeAccelValues(stacks)
armLengths = getArmLengths(stacks)

# plt.figure(0)
# plt.plot(newTimes, armLengths[0])
# plt.plot(newTimes, armLengths[1])
# plt.plot(newTimes, armLengths[2])


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

plt.figure(2)
plt.plot(newTimes, armDif[0], label = 'arm 1 displacement')
plt.plot(newTimes, armDif[1], label = 'arm 2 displacement')
plt.plot(newTimes, armDif[2], label = 'arm 3 displacement')
plt.grid()
plt.xlabel('time (years)')
plt.ylabel('displacement (meters)')
plt.title('Perturbation of LISA Arm Lengths due to Ceres\nat 0 Degree Initial Phase Angle')#\nwith 0 constellation phase angle')
plt.legend()
plt.ylim([-20, 20])

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
