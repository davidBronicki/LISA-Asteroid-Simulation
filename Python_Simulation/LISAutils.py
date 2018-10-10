'''This script revolves around the orbit class. All external variables and
methods are implemented to support the orbit class.
This class holds orbital parameters and allows for simple propogation of orbits.
Note:
When setting mean anomaly or time, use setMean, setTime, or iterateTime. These
methods will automatically change the other anomalies accordingly.
Variables are:
f = true anomaly
psi = eccentric anomaly
eta = mean anomaly
e = eccentricity
i = inclination of orbital plane from ecliptic plane
omega = angle of ascending node from x axis
theta = angle of periapsis from ascending node
p = semi latus rectum
a = semi major axis
t = time
mass = mass
Also:
t0 and eta0 are initial values
co, so are sine and cosine of omega
ci, si are sine and cosine of inclination
these are made to save on computation power (not really important)'''
import numpy as np
from math import *

G = 6.67408 * 10 ** -11
solarMass = 1.989 * 10 ** 30
earthMass = 5.972 * 10 ** 24
astroUnit = 149597870700

mu = G * solarMass

def yearToEpochMJD(inputYear):
    return 365 * (inputYear - 2030) + 62502

'''The following group of functions are for converting from
one type of anomaly to the other. While this is only implemented
for use in the orbit class, they are defined in global scope
so they may be used for general purpose.'''

def TrueToEccAnom(f, ecc):
    temp = sqrt((1 - ecc) / (1 + ecc)) * tan(f / 2)
    return (2 * np.arctan(temp)) % (2 * pi)

def EccToTrueAnom(psi, ecc):
    temp = sqrt((1 + ecc) / (1 - ecc)) * tan(psi / 2)
    return (2 * np.arctan(temp)) % (2 * pi)

def EccToMeanAnom(psi, ecc):
    return psi - ecc * sin(psi)

def DifOfMean(eta, psi, ecc):
    return eta - psi + ecc * sin(psi)

def MeanToEccAnom(eta, ecc):
    eta = eta % (2 * pi)
    width = 2 * pi
    lower = 0
    upper = 2 * pi
    midpoint = pi
    while (width > 10 ** -14):
        if DifOfMean(eta, midpoint, ecc) * DifOfMean(eta, lower, ecc) > 0:
            lower = midpoint
        else:
            upper = midpoint
        width = upper - lower
        midpoint = (upper + lower) / 2
    return midpoint

class orbit:
    '''There are four optional keywords to initialize this class.

    anomType can be:
    trueAnom: means the first argument represents true anomaly
    meanAnom: first argument is mean anomaly
    eccAnom: first argument is eccentric anomaly

    dimType can be:
    semiMajor: means the 'dimension' argument represents semi major axis
    semiLatus: 'dimension' represents semi latus rectum

    time is the initial time. My thinking is that this may be useful
    for creating a new orbit on the fly after the rest of the system
    has already started. May need revision to mesh with the standard
    time variable used in the field (the epoch).

    mass is self explanatory. Not needed for a satalite but is for a planet.'''
    
    def __init__(self, anomaly, ecc, incline, angAscend, angPeri, dimension,
                anomType = 'trueAnom', dimType = 'semiMajor', simTime = 0,
                paramTime = 0, mass = 1, name = 'defaultName'):
        self.e = ecc
        if anomType == 'meanAnom':
            self.eta = anomaly % (2 * pi)
            self.eta0 = self.eta
            self.psi = MeanToEccAnom(self.eta, ecc)
            self.f = EccToTrueAnom(self.psi, ecc)
        elif anomType == 'eccAnom':
            self.psi = anomaly % (2 * pi)
            self.eta = EccToMeanAnom(self.psi, ecc)
            self.eta0 = self.eta
            self.f = EccToTrueAnom(self.psi, ecc)
        else: #type = true anomaly (default)
            self.f = anomaly % (2 * pi)
            self.psi = TrueToEccAnom(self.f, ecc)
            self.eta = EccToMeanAnom(self.psi, ecc)
            self.eta0 = self.eta
        
        self.i = incline
        self.omega = angAscend
        self.theta = angPeri
        
        if dimType == 'semiLatus':
            self.p = dimension
            self.a = self.p / (1 - ecc * ecc)
        else: #type = semi major axis (default)
            self.a = dimension
            self.p = self.a * (1 - ecc * ecc)
        
        self.meanAngMotion = sqrt(mu / self.a ** 3)
        self.t = 0
        self.eta = (self.eta0 + (simTime - paramTime)
            * self.meanAngMotion * 24 * 3600) % (2 * np.pi)
        self.eta0 = self.eta
        self.mass = mass
        self.name = name

        #commonly used variables (saves on processing power)
        self.co = cos(self.omega)
        self.so = sin(self.omega)
        self.ci = cos(self.i)
        self.si = sin(self.i)
        if self.si == 0:
            self.si = 10 ** -10
        
        self.makeR()

    def __str__(self):
        params = 'Orbital Parameters:\n' + \
            '\tSemi-Major Axis: ' + str(degrees(self.a)) + '\n' + \
            '\tEccentricity: ' + str(degrees(self.e)) + '\n' + \
            '\tInclination: ' + str(degrees(self.i)) + '\n' + \
            '\tLongitude of Ascending Node: ' + str(degrees(self.omega)) + '\n' + \
            '\tArgument of Perihelion: ' + str(degrees(self.theta)) + '\n' + \
            '\tCurrent True Anomaly: ' + str(degrees(self.f))
        name = 'Name: ' + self.name
        mass = 'Mass: ' + str(self.mass)
        return name + '\n' + mass + '\n' + params

    def makeR(self):
        self.r = self.p / (1 + self.e * cos(self.f))

    def iterateTime(self, dt):#set time and iterate time call setMean
        self.setTime(self.t + dt)

    def setMean(self, newEta):#also calculates eccentric and true anomalies and new time
        self.eta = newEta % (2 * pi)
        self.psi = MeanToEccAnom(self.eta, self.e)
        self.f = EccToTrueAnom(self.psi, self.e)
        self.makeR()

    def setTime(self, newT):
        self.t = newT
        self.setMean(newT * self.meanAngMotion + self.eta0)

    def x(self):
        return self.r * (self.co * cos(self.theta + self.f)
            - self.ci * self.so * sin(self.theta + self.f))
    def y(self):
        return self.r * (self.so * cos(self.theta + self.f)
            + self.ci * self.co * sin(self.theta + self.f))
    def z(self):
        return self.r * (self.si * sin(self.theta + self.f))
    def pos(self):
        return np.array([self.x(), self.y(), self.z()])

    def vx(self):
        return -sqrt(mu/self.p)*(self.co*(sin(self.theta+self.f)+self.e*sin(self.theta))
                                 +self.ci*self.so*(cos(self.theta+self.f)+self.e*cos(self.theta)))
    def vy(self):
        return -sqrt(mu/self.p)*(self.so*(sin(self.theta+self.f)+self.e*sin(self.theta))
                                 -self.ci*self.co*(cos(self.theta+self.f)+self.e*cos(self.theta)))
    def vz(self):
        return sqrt(mu/self.p)*self.si*(cos(self.theta+self.f)+self.e*cos(self.theta))
    def vel(self):
        return np.array([self.vx(), self.vy(), self.vz()])

    def forceFrom(self, otherBody):
        distance = dist(self, otherBody)
        gamma = G * self.mass * otherBody.mass
        disp = displacement(self, otherBody)
        return -disp * gamma / distance**3
    def forceOn(self, otherBody):
        distance = dist(self, otherBody)
        gamma = G * self.mass * otherBody.mass
        disp = displacement(self, otherBody)
        return disp * gamma / distance**3

    '''I'm not sure the perturbation section is working, but I created the methods
    based on the perturbation section of the book packet you gave me.

    There is a method for perturbing from another body (orbit class) which creates
    a force and passes it to the more general perturbation method. In the general perturbation
    method, the force is first converted to cylindrical coordinates, then the time dependent
    perturbation equations are used to change all of the orbital parameters. Keep in mind that
    the time dependent equations do not use the first order approximation that the true anomaly
    dependant equations use. I don't know how big of an effect that has.

    Aside from that, beginPerturb should be called first since it makes what is effectively the
    coordinate transform matrix that gets used to convert the force into cylindrical.
    Once everything is perturbed, endPerturb must be called since it sets variables that
    are used in the position and velocity methods.'''
    def beginPerturb(self):#sets up variables for perturbation
        ctotal = cos(self.f + self.theta)
        stotal = sin(self.f + self.theta)

        #r,phi,z unit vectors in cartesian coordinates
        
        self.rx = self.co * ctotal - self.ci * self.so * stotal
        self.ry = self.so * ctotal + self.ci * self.co * stotal
        self.rz = self.si * stotal

        self.phix = -self.co * stotal - self.ci * self.so * ctotal
        self.phiy = -self.so * stotal + self.ci * self.co * ctotal
        self.phiz = self.si * ctotal

        self.zx = self.si * self.so
        self.zy = -self.si * self.co
        self.zz = self.ci
    
    def perturb(self, fx, fy, fz, dt):
        fR = fx * self.rx + fy * self.ry + fz * self.rz
        fPhi = fx * self.phix + fy * self.phiy + fz * self.phiz
        fZ = fx * self.zx + fy * self.zy + fz * self.zz
        
        commonSqrt = sqrt(self.p / mu)
        commonDen = 1 / (1 + self.e * cos(self.f))
        
        dp = dt * 2 * self.p * commonSqrt * fPhi * commonDen
        de = dt * commonSqrt * (sin(self.f) * fR
                + (2 * cos(self.f) + self.e * (1 + cos(self.f) ** 2)) * fPhi * commonDen)
        dIncline = dt * commonSqrt * cos(self.f + self.theta) * fZ * commonDen
        dOmega = dt * commonSqrt * sin(self.f + self.theta) * fZ * commonDen / self.si
        dTheta = dt * commonSqrt * (-cos(self.f) * fR
                + (2 + self.e * cos(self.f)) * sin(self.f) * fPhi * commonDen
                - self.e * (self.ci / self.si) * sin(self.theta + self.f) * fZ * commonDen) / self.e
        df = -dTheta-self.ci*dOmega
        
        self.p += dp
        self.e += de
        self.i += dIncline
        self.omega += dOmega
        self.theta += dTheta
        self.f += df

    def perturbFrom(self, other, dt):#sets up forces from another object and applies them
        distx = self.x() - other.x()
        disty = self.y() - other.y()
        distz = self.z() - other.z()
        distance = sqrt(distx * distx + disty * disty + distz * distz)
        temp = -G * other.mass / (distance ** 3)
        self.perturb(distx * temp, disty * temp, distz * temp, dt)

    def endPerturb(self):#reassigns commonly used variables
        self.a = self.p / (1 - self.e*self.e)
        self.co = cos(self.omega)
        self.so = sin(self.omega)
        self.ci = cos(self.i)
        self.si = sin(self.i)
        if self.si == 0:
            self.si = 10 ** -10

# def stateSpaceOrbit(pos, vel, time = 0, mass = 1):
#     pass

def displacement(fromOrbit, toOrbit):
    return toOrbit.pos() - fromOrbit.pos()

def distSqr(orbit1, orbit2):
    sqr = (orbit1.x() - orbit2.x()) ** 2
    sqr += (orbit1.y() - orbit2.y()) ** 2
    sqr += (orbit1.z() - orbit2.z()) ** 2
    return sqr

def dist(orbit1, orbit2):
    return sqrt(distSqr(orbit1, orbit2))


class LISA:
    armLength = 2.5*10**9
    def __init__(self, initialOrientiationAngle = 0, initialLongitude = 0):
        alpha = LISA.armLength / (2 * astroUnit)
        ecc = sqrt(1 + (2 / sqrt(3)) * alpha + (4 / 3) * alpha ** 2) - 1
        incline = np.arctan(alpha / (1 + alpha / sqrt(3)))

        orbitalOffsetAngle = 2 * pi / 3#each orbit is a 120 degree rotation of the others

        self.sats = []

        for i in range(3):
            self.sats.append(orbit(
                pi / 2 + initialOrientiationAngle + i * orbitalOffsetAngle + initialLongitude,             #mean anomaly
                ecc,                                                                    #ecc
                incline,                                                                #incline
                - initialOrientiationAngle - i * orbitalOffsetAngle,   #long. of asc. node
                -pi / 2,                                                                #arg. of Peri.
                astroUnit,                                                              #semi major axis
                anomType = 'meanAnom', name = 'LISA Craft ' + str(i)))

    def sat(self, n):
        return self.sats[n]

    def setTime(self, time):
        for sat in self.sats:
            sat.setTime(time)

    def pos(self):
        positions = []
        for sat in self.sats:
            positions.append(sat.pos())
        return positions

    def vel(self):
        velocities = []
        for sat in self.sats:
            velocities.append(sat.vel())
        return velocities

    def getAccelFrom(self, otherBody):
        output = []
        for sat in self.sats:
            output.append(sat.forceFrom(otherBody))
        return np.array(output)

    def armLengths(self):
        positions = self.pos()
        lengths = []
        for i in range(3):
            for j in range(i+1, 3):
                r = positions[i] - positions[j]
                lengths.append(np.linalg.norm(r))
        return lengths

    def distancesToObject(self, otherObject):
        otherPosition = otherObject.pos()
        lengths = []
        for sat in self.sats:
            r = sat.pos() - otherPosition
            lengths.append(np.linalg.norm(r))
        return lengths

    def beginPerturb(self):
        for sat in self.sats:
            sat.beginPerturb()

    def perturb(self, force, dt):
        for sat in self.sats:
            sat.perturb(force[0], force[1], force[2], dt)

    def pertrubFrom(self, otherObject, dt):
        for sat in self.sats:
            sat.perturbFrom(otherObject, dt)

    def endPerturb(self):
        for sat in self.sats:
            sat.endPerturb()

    def units(self):
        locations = self.pos()
        unit12 = (locations[1] - locations[0])/np.linalg.norm(locations[1]-locations[0])
        unit13 = (locations[2] - locations[0])/np.linalg.norm(locations[2]-locations[0])
        unit23 = (locations[2] - locations[1])/np.linalg.norm(locations[2]-locations[1])
        return [unit12, unit13, unit23]

    def psuedoPerturbForceFromList(self, objects):
        a = []
        locations = self.pos()
        for i in range(3):
            a.append(np.array([0,0,0], dtype = 'float'))
        for otherObject in objects:
            otherPos = otherObject.pos()
            for i in range(3):
                displacement = locations[i] - otherPos
                temp = -G * otherObject.mass / np.linalg.norm(displacement)**3
                a[i] += displacement * temp

        units = self.units()
        
        discrepency12 = np.dot(a[1]-a[0], units[0])
        discrepency13 = np.dot(a[2]-a[0], units[1])
        discrepency23 = np.dot(a[2]-a[1], units[2])

        return [discrepency12, discrepency13, discrepency23]

    def createForceDiscrepancyList(self, objects):
        locations = self.pos()
        out12 = []
        out13 = []
        out23 = []
        r12 = (locations[1] - locations[0])/np.linalg.norm(locations[1] - locations[0])
        r13 = (locations[2] - locations[0])/np.linalg.norm(locations[2] - locations[0])
        r23 = (locations[2] - locations[1])/np.linalg.norm(locations[2] - locations[1])
        for thing in objects:
            pos = thing.pos()
            R1 = locations[0] - pos
            R2 = locations[1] - pos
            R3 = locations[2] - pos
            out12.append(np.dot(R2/np.linalg.norm(R2)**3 - R1/np.linalg.norm(R1)**3, r12))
            out13.append(np.dot(R3/np.linalg.norm(R3)**3 - R1/np.linalg.norm(R1)**3, r13))
            out23.append(np.dot(R3/np.linalg.norm(R3)**3 - R2/np.linalg.norm(R2)**3, r23))
        return np.array([out12,out13, out23])

    def areaChangeRatio(self, discrepencies):
        r1 = self.sats[0].pos()
        r2 = self.sats[1].pos()
        r3 = self.sats[2].pos()
        a = np.linalg.norm(r2 - r1)
        b = np.linalg.norm(r3 - r1)
        c = np.linalg.norm(r3 - r2)
        s = (a + b + c) / 2
        #A = sqrt(s*(s-a)*(s-b)*(s-c))
        da = discrepencies[0]
        db = discrepencies[1]
        dc = discrepencies[2]
        #returns dA / A
        return (1 / 4) * (da * (1/s - 1/(s-a) + 1/(s-b) + 1/(s-c)) +
                          db * (1/s + 1/(s-a) - 1/(s-b) + 1/(s-c)) +
                          dc * (1/s + 1/(s-a) + 1/(s-b) - 1/(s-c)))

