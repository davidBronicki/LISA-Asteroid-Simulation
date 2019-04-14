import my_package.util as tools
import numpy as np
import LISAutils

gaussianValues = []
gaussCoef = 1/np.sqrt(2*np.pi)
currentGaussIndex = 0

def gaussDist(x):
	return gaussCoef * np.exp(-x**2/2)

def initGaussianValues(initNumber):
	gaussianValues = tools.metropolisAlgorithm(gaussDist, initNumber, 0, dryRunCount = 5)

def getGaussianValue(mean, sigma):
	x = gaussianValues[currentGaussIndex]
	currentGaussIndex = (127 + currentGaussIndex) % len(gaussianValues)
	return x * sigma + mean

def interpSpectralType(inputString):
	if (inputString == "C" or inputString == "B"
		or inputString == "F" or inputString == "F"
		or inputString == "G" or inputString == "D"
		or inputString == "T"):
		return "Carb"
	elif (inputString == "S" or inputString == "V"
		or inputString == "A" or inputString == "R"):
		return "Stony"
	elif (inputString == "M"):
		return "Metal"
	else: return "Other"

class orbitalDataObject:
	def __init__(self, dataStringList):
		self.name = dataStringList[0]
		self.epoch = float(dataStringList[1])
		self.orbitalParams = []
		for i in range(2, 8):
			self.orbitalParams.append(float(dataStringList[i]))
		self.sigmaOrbitalParams = []
		for i in range(8, 14):
			self.sigmaOrbitalParams.append(float(dataStringList[i]))
		try:
			self.mass = float(dataStringList[14])
		except:
			self.mass = 0
		self.absMag = float(dataStringList[15])
		try:
			self.sigmaAbsMag = float(dataStringList[16])
		except:
			self.sigmaAbsMag = 0
		try:
			self.albedo = float(dataStringList[17])
		except:
			self.albedo = 0
		self.specType = interpSpectralType(dataStringList[18])

	def findMass(self, absMag, type, albedo):


	def getOrbitalObject(self):
		orbitalParams = []
		for indexPair in zip(self.orbitalParams, self.sigmaOrbitalParams):
			orbitalParams.append(getGaussianValue(indexPair[0], indexPair[1]))
		if (self.sigmaAbsMag == 0):
			absMag = getGaussianValue(self.absMag, .2)
		else: absMag = getGaussianValue(self.absMag, self.sigmaAbsMag)
		if (self.mass == 0):
			mass = findMass(absMag, self.type, self.albedo)
		else:
			mass = self.mass / 6.67e-2
		return LISAutils.Orbit(orbitalParams[5], orbitalParams[1], orbitalParams[2],
			orbitalParams[3], orbitalParams[4], orbitalParams[0], 'meanAnom',
			)

def extractData(fileLocation, numberToUpload = 0):
	file = open(fileLocation, 'r')
	n = 0
	outputData = []
	for line in file:
		if numberToUpload == 0 or n < numberToUpload:
			try:
				temp = orbitalDataObject(line)
				outputData.append(temp)
				n += 1
			except:
				pass
		else:
			break
	return outputData[1:]
