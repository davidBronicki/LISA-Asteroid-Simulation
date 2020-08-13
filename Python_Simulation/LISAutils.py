'''
This script revolves around the orbit class. All external variables and
methods are implemented to support the orbit class.
This class holds orbital parameters and allows for simple time evolution of orbits.
'''

import numpy as np
from math import sqrt, sin, cos, tan, pi

G = 6.67408 * 10 ** -11
solarMass = 1.989 * 10 ** 30
earthMass = 5.972 * 10 ** 24
astroUnit = 149597870700

mu = G * solarMass

def yearToEpochMJD(inputYear):
    return 365 * (inputYear - 2030) + 62502

'''
The following group of functions are for converting from
one type of anomaly to another. While this is only implemented
for use in the orbit class, they are defined in global scope
so they may be used for general purpose.
'''

def TrueToEccAnom(f, ecc):
    temp = sqrt((1 - ecc) / (1 + ecc)) * tan(f / 2)
    return (2 * np.arctan(temp)) % (2 * pi)

def EccToTrueAnom(psi, ecc):
    temp = sqrt((1 + ecc) / (1 - ecc)) * tan(psi / 2)
    return (2 * np.arctan(temp)) % (2 * pi)

def EccToMeanAnom(psi, ecc):
    return psi - ecc * sin(psi)

def MeanToEccAnom(eta, ecc):
    def DifOfMean(eta, psi, ecc):
        return eta - psi + ecc * sin(psi)
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
    '''
    Note:
    When setting mean anomaly or time, use setMean, setTime, or iterateTime. These
    methods will automatically change the other anomalies accordingly.

    IMPORTANT:
    There is a bug when mixing 'setTime' and 'iterateTime' in code. They seem
    to work independently, but don't play well together for some reason.
    So avoid mixing and matching when you actually use this.

    Variables are:
        class variable name (initializer name) = description of the variable

        f (anomaly) = true anomaly
        psi (anomaly) = eccentric anomaly
        eta (anomaly) = mean anomaly
        e (ecc) = eccentricity
        i (incline) = inclination of orbital plane from ecliptic plane
        RAAN (angAscend) = Right Assention of the Assending Node (capital omega in literature)
        omega (angPeri) = angle to the periapsis from ascending node (lower case omega in literature)
        p (dimension) = semi latus rectum
        a (dimension) = semi major axis
        t (simTime and paramTime, see below) = time
        mass = mass
    Also:
        t0 and eta0 are initial values
        co, so are sine and cosine of RAAN ('o' is for capital Omega)
        ci, si are sine and cosine of inclination
            These are made to save on computation power (not important for the physics)
        and various variables are defined for perturbation purposes (see perturbation functions)


    There are six optional keywords to initialize this class.

    anomType can be:
        trueAnom: means the 'anomaly' argument represents true anomaly
        meanAnom: 'anomaly' represents mean anomaly
        eccAnom: 'anomaly' represents eccentric anomaly

    dimType can be:
        semiMajor: means the 'dimension' argument represents semi major axis
        semiLatus: 'dimension' represents semi latus rectum

    simTime is the MJD time that the simulation starts

    paramTime is the MJD epoch that is given with the orbital elements

    mass and name are self explanatory.

    Different times used:
        t will always equal zero at the beginning of the simulation;
        this is the number of seconds since the beginning of the simulation.

        simTime is an initialization parameter marking the MJD time
        (number of days since midnight of November 17th, 1858)
        at which the simulation starts. If you want a physically accurate
        simulation of a real world object, you should either find a cite
        to generate the MJD of when your simulation starts or
        you should have a piece of code to generate it.

        paramTime is an initialization parameter marking the MJD time
        of the orbital parameters provided. This should be given with
        the orbital parameters, often labelled 'epoch'. For instance,
        JPL gives this under "Epoch (MJD)" (there are multiple ways of
        giving the Epoch, so the method is given in parentheses).

        During initialization, the orbital parameters are NON-PERTURBATIVELY
        time evolved from "paramTime" to "simTime", and then these are locked
        in. (This amounts to only evolving the anomalies, since all
        other parameters are constant when there are no perturbations.)
        Since this is done non-perturbatively, the simulation will
        lose fidelity if "simTime" is too far in the future
        or past from "paramTime".
    '''
    
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
        self.RAAN = angAscend
        self.omega = angPeri
        
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
        self.co = cos(self.RAAN)
        self.so = sin(self.RAAN)
        self.ci = cos(self.i)
        self.si = sin(self.i)
        if self.si == 0:
            self.si = 10 ** -10
        
        self.psi = MeanToEccAnom(self.eta, self.e)
        self.f = EccToTrueAnom(self.psi, self.e)
        self.makeR()

    def __str__(self):
        params = 'Orbital Parameters:\n' + \
            '\tSemi-Major Axis: ' + str(degrees(self.a)) + '\n' + \
            '\tEccentricity: ' + str(degrees(self.e)) + '\n' + \
            '\tInclination: ' + str(degrees(self.i)) + '\n' + \
            '\tLongitude of Ascending Node: ' + str(degrees(self.RAAN)) + '\n' + \
            '\tArgument of Perihelion: ' + str(degrees(self.omega)) + '\n' + \
            '\tCurrent True Anomaly: ' + str(degrees(self.f))
        name = 'Name: ' + self.name
        mass = 'Mass: ' + str(self.mass)
        return name + '\n' + mass + '\n' + params

    def makeR(self):
        self.r = self.p / (1 + self.e * cos(self.f))

    def iterateTime(self, dt):#step forward through time
        self.t+= dt
        self.setMean(self.eta + dt * self.meanAngMotion)
        # self.setTime(self.t + dt)

    def setMean(self, newEta):#also calculates eccentric and true anomalies and new time
        self.eta = newEta % (2 * pi)
        self.psi = MeanToEccAnom(self.eta, self.e)
        self.f = EccToTrueAnom(self.psi, self.e)
        self.makeR()

    def setTime(self, newT):
        self.t = newT
        self.setMean(newT * self.meanAngMotion + self.eta0)

    def x(self):
        return self.r * (self.co * cos(self.omega + self.f)
            - self.ci * self.so * sin(self.omega + self.f))
    def y(self):
        return self.r * (self.so * cos(self.omega + self.f)
            + self.ci * self.co * sin(self.omega + self.f))
    def z(self):
        return self.r * (self.si * sin(self.omega + self.f))
    def pos(self):
        return np.array([self.x(), self.y(), self.z()])

    def vx(self):
        return -sqrt(mu/self.p)*(self.co*(sin(self.omega+self.f)+self.e*sin(self.omega))
            +self.ci*self.so*(cos(self.omega+self.f)+self.e*cos(self.omega)))
    def vy(self):
        return -sqrt(mu/self.p)*(self.so*(sin(self.omega+self.f)+self.e*sin(self.omega))
            -self.ci*self.co*(cos(self.omega+self.f)+self.e*cos(self.omega)))
    def vz(self):
        return sqrt(mu/self.p)*self.si*(cos(self.omega+self.f)+self.e*cos(self.omega))
    def vel(self):
        return np.array([self.vx(), self.vy(), self.vz()])

    '''
    These calculate the force (not force per unit mass) which another
    body exerts on this body, and which this body exerts on another
    body, respectively. (They, of course, only differ by a minus sign.)
    '''
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

    '''
    beginPerturb should be called first, and once everything is perturbed,
    endPerturb must be called to set miscellaneous variables
    not set in the actual perturbation code

    The perturb function accepts the force in x y z coordinates, and the time interval used.

    The perturbFrom function accepts another orbit class (which has an optional mass argument)
    and calculates the acceleration felt from this object, then passes this to
    the perturb function.
    '''
    def beginPerturb(self):#sets up variables for perturbation
        ctotal = cos(self.f + self.omega)
        stotal = sin(self.f + self.omega)

        #building orbit based unit vectors in xyz coordinates
        
        #r is the typical radial unit vector
        self.rx = self.co * ctotal - self.ci * self.so * stotal
        self.ry = self.so * ctotal + self.ci * self.co * stotal
        self.rz = self.si * stotal

        #phi is the unit vector orthogonal to r in the plane of orbit
        #pointing in the prograde direction
        self.phix = -self.co * stotal - self.ci * self.so * ctotal
        self.phiy = -self.so * stotal + self.ci * self.co * ctotal
        self.phiz = self.si * ctotal

        #z is the unit vector orthogonal to the plane of orbit with z=r cross phi
        #(to ensure a right handed coordinate system)
        self.zx = self.si * self.so
        self.zy = -self.si * self.co
        self.zz = self.ci
    
    def perturb(self, fx, fy, fz, dt):
        #Remember to call begin perturb before this and end perturb after.
        #If there are several perturbing forces, you only need to call
        #begin once and end once.

        #essentially a coordinate change for the force.
        #change from x,y,z unit vectors to orbit based unit vectors (see above function)
        fR = fx * self.rx + fy * self.ry + fz * self.rz
        fPhi = fx * self.phix + fy * self.phiy + fz * self.phiz
        fZ = fx * self.zx + fy * self.zy + fz * self.zz

        #basically black box equations for how the orbital elements
        #change due to an acceleration over a time dt (aka a small impulse)
        adot = 2*sqrt(self.a**3/mu/(1-self.e**2))*(self.e * fR*sin(self.f)
            + fPhi *(1+self.e * cos(self.f)))
        edot = sqrt(self.a*(1-self.e**2)/mu)*(fR*sin(self.f)+fPhi
            *(cos(self.f)+cos(self.psi)))
        idot = sqrt(self.a*(1-self.e**2)/mu)*fZ*cos(self.f+self.omega)/(1+self.e*cos(self.f))
        RAANdot = sqrt(self.a*(1-self.e**2)/mu)*fZ*sin(self.f+self.omega)\
            /sin(self.i)/(1+self.e*cos(self.f))
        omegadot = -RAANdot*cos(self.i)+sqrt(self.a*(1-self.e**2)/self.e**2/mu)\
            *(-fR * cos(self.f)+fPhi * (2+self.e*cos(self.f))*sin(self.f)
                /(1+self.e*cos(self.f)))
        fDot = -omegadot-self.ci*RAANdot

        #apply the perturbations
        self.a += dt*adot
        self.e += dt*edot
        self.i += dt*idot
        self.RAAN += dt*RAANdot
        self.omega += dt*omegadot
        self.f += dt*fDot


    def perturbFrom(self, other, dt):
        #sets up forces from another object and sends it to the perturb function
        distx = self.x() - other.x()
        disty = self.y() - other.y()
        distz = self.z() - other.z()
        distance = sqrt(distx * distx + disty * disty + distz * distz)
        temp = -G * other.mass / (distance ** 3)
        self.perturb(distx * temp, disty * temp, distz * temp, dt)

    def endPerturb(self):
        #reassigns commonly used variables
        self.p = self.a * (1-self.e**2)
        self.co = cos(self.RAAN)
        self.so = sin(self.RAAN)
        self.ci = cos(self.i)
        self.si = sin(self.i)
        if self.si == 0:
            self.si = 10e-10

        oldMeanAngMotion = self.meanAngMotion
        self.meanAngMotion = sqrt(mu / self.a ** 3)
        deltaMeanAngMotion = self.meanAngMotion - oldMeanAngMotion
        self.eta0 -= deltaMeanAngMotion * self.t

        oldEta = self.eta

        self.f = self.f % (2 * pi)
        self.psi = TrueToEccAnom(self.f, self.e)
        self.eta = EccToMeanAnom(self.psi, self.e)

        deltaEta = self.eta - oldEta
        self.eta0 += deltaEta

        self.makeR()

#These functions should be pretty self explanatory:
#quick ways to get distances between orbiting objects.

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
    '''
    This is a specialized class for simulating the orbits of the LISA crafts.

    initialOrientiationAngle is the 'rotation' of the LISA triangle
    initialLongitude is the initial right ascension of the 'center' of the LISA triangle

    During initialization, three "Orbit" objects are made, one for each LISA craft.
    The exact conditions are set according to the initialization variables
    so that at t = 0 (the start of the simulation), the crafts satisfy
    the desired conditions.

    Bear in mind that "initialOrientationAngle" differs from the standard
    by pi/3 radians. So if you would like to initialize with an angle of theta,
    then pass (theta - pi/3) to the initializer. (Perhaps someone can fix this
    at some point, I'm just too lazy to do it.)
    '''
    armLength = 2.5e9
    def __init__(self, initialOrientiationAngle = 0, initialLongitude = 0, flipLISA = False):
        alpha = LISA.armLength / (2 * astroUnit)
        ecc = sqrt(1 + (2 / sqrt(3)) * alpha + (4 / 3) * alpha ** 2) - 1
        incline = np.arctan(alpha / (1 + alpha / sqrt(3)))

        orbitalOffsetAngle = 2 * pi / 3#each orbit is a 120 degree rotation of the others

        self.sats = []

        if flipLISA:
            for i in range(3):
                self.sats.append(orbit(
                    pi / 2 + initialOrientiationAngle\
                        + i * orbitalOffsetAngle + initialLongitude,        #mean anomaly
                    ecc,                                                    #ecc
                    incline,                                                #incline
                    -initialOrientiationAngle - i * orbitalOffsetAngle + pi,#long. of asc. node
                    pi / 2,                                                 #arg. of Peri.
                    astroUnit,                                              #semi major axis
                    anomType = 'meanAnom', name = 'LISA Craft ' + str(i)))
        else:
            for i in range(3):
                self.sats.append(orbit(
                    pi / 2 + initialOrientiationAngle\
                        + i * orbitalOffsetAngle + initialLongitude,        #mean anomaly
                    ecc,                                                    #ecc
                    incline,                                                #incline
                    - initialOrientiationAngle - i * orbitalOffsetAngle,    #long. of asc. node
                    -pi / 2,                                                #arg. of Peri.
                    astroUnit,                                              #semi major axis
                    anomType = 'meanAnom', name = 'LISA Craft ' + str(i)))
                tempPos = self.sats[-1].pos()
                azimuth = np.arctan2(tempPos[1],tempPos[0])
                polar = np.arctan(tempPos[2]/np.sqrt(tempPos[1]**2 + tempPos[0]**2))

        tempPos = self.pos()
        tempPosAvg = sum(tempPos)/3
        normal = np.cross(tempPos[1] - tempPos[0], tempPos[1] - tempPos[2])
        normal = normal / np.linalg.norm(normal)
        print(tempPosAvg)
        print(normal)

    def sat(self, n):
        return self.sats[n]

    def iterateTime(self, dt):#step forward through time
        for sat in self.sats:
            sat.iterateTime(dt)

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

    def perturbFrom(self, otherObject, dt):
        for sat in self.sats:
            sat.perturbFrom(otherObject, dt)

    def endPerturb(self):
        for sat in self.sats:
            sat.endPerturb()

    '''
    The following functions may be useful, but it has been a long time
    since I used them and I have forgotten what they do.
    Hence they are commented out; I cannot guarantee they work properly
    nor can I explain what they do.
    '''

    # def units(self):
    #     locations = self.pos()
    #     unit12 = (locations[1] - locations[0])/np.linalg.norm(locations[1]-locations[0])
    #     unit13 = (locations[2] - locations[0])/np.linalg.norm(locations[2]-locations[0])
    #     unit23 = (locations[2] - locations[1])/np.linalg.norm(locations[2]-locations[1])
    #     return [unit12, unit13, unit23]

    # def psuedoPerturbForceFromList(self, objects):
    #     a = []
    #     locations = self.pos()
    #     for i in range(3):
    #         a.append(np.array([0,0,0], dtype = 'float'))
    #     for otherObject in objects:
    #         otherPos = otherObject.pos()
    #         for i in range(3):
    #             displacement = locations[i] - otherPos
    #             temp = -G * otherObject.mass / np.linalg.norm(displacement)**3
    #             a[i] += displacement * temp

    #     units = self.units()
        
    #     discrepency12 = np.dot(a[1]-a[0], units[0])
    #     discrepency13 = np.dot(a[2]-a[0], units[1])
    #     discrepency23 = np.dot(a[2]-a[1], units[2])

    #     return [discrepency12, discrepency13, discrepency23]

    # def createForceDiscrepancyList(self, objects):
    #     locations = self.pos()
    #     out12 = []
    #     out13 = []
    #     out23 = []
    #     r12 = (locations[1] - locations[0])/np.linalg.norm(locations[1] - locations[0])
    #     r13 = (locations[2] - locations[0])/np.linalg.norm(locations[2] - locations[0])
    #     r23 = (locations[2] - locations[1])/np.linalg.norm(locations[2] - locations[1])
    #     for thing in objects:
    #         pos = thing.pos()
    #         R1 = locations[0] - pos
    #         R2 = locations[1] - pos
    #         R3 = locations[2] - pos
    #         out12.append(np.dot(R2/np.linalg.norm(R2)**3 - R1/np.linalg.norm(R1)**3, r12))
    #         out13.append(np.dot(R3/np.linalg.norm(R3)**3 - R1/np.linalg.norm(R1)**3, r13))
    #         out23.append(np.dot(R3/np.linalg.norm(R3)**3 - R2/np.linalg.norm(R2)**3, r23))
    #     return np.array([out12,out13, out23])

    # def areaChangeRatio(self, discrepencies):
    #     r1 = self.sats[0].pos()
    #     r2 = self.sats[1].pos()
    #     r3 = self.sats[2].pos()
    #     a = np.linalg.norm(r2 - r1)
    #     b = np.linalg.norm(r3 - r1)
    #     c = np.linalg.norm(r3 - r2)
    #     s = (a + b + c) / 2
    #     #A = sqrt(s*(s-a)*(s-b)*(s-c))
    #     da = discrepencies[0]
    #     db = discrepencies[1]
    #     dc = discrepencies[2]
    #     #returns dA / A
    #     return (1 / 4) * (da * (1/s - 1/(s-a) + 1/(s-b) + 1/(s-c)) +
    #                       db * (1/s + 1/(s-a) - 1/(s-b) + 1/(s-c)) +
    #                       dc * (1/s + 1/(s-a) + 1/(s-b) - 1/(s-c)))

