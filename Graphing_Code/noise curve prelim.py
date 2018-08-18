import numpy as np
import matplotlib.pyplot as plt

file = open('lisa_noise_curve.txt')

data = []

for line in file:
	temp = line.split(' ')
	try:
		tempData = []
		for item in temp:
			tempData.append(float(item))
		data.append(tempData)
	except:
		pass

data = list(zip(*data))

freq = np.array(list(data[0]))
strain = np.array(list(data[1]))



plt.plot(np.log10(freq), np.log10(strain))
plt.grid()
plt.show()

