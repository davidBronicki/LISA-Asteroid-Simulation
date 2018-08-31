import numpy as np

minMaxValues = [29.58,30.19,30.7,31.31,31.85,32.39,33.02,33.55,33.81,34.11]
closeApproachValues = [29.36,30.61,31.89,33.22]
# minMaxValues = [29.58,30.19,30.7,31.31,31.85,32.39,33.02]
# closeApproachValues = [29.36,30.61,31.89]

output = []
for approach in closeApproachValues:
	temp = []
	for item in minMaxValues:
		temp.append(item - approach)
	output.append(temp)
output = np.array(output)
print(output)
