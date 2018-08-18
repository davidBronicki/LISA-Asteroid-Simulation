import custum.usefulTools as tools
import numpy as np
import matplotlib.pyplot as plt

root = R"../../reference_file/data_dumps/"
postfix = '.csv'

# extraBit = 'Relative Acceleration '
extraBit = 'Absolute Acceleration '
# extraBit = 'Unit Vector '

# label = '-10k Asteroids_-30k Points of Time Sampled'
# label = '-Jupiter_-10k Points of Time Sampled -10 years'
# label = '-Ceres_-10k Points of Time Sampled'
# label = '-     4 Vesta_-10k Points of Time Sampled'
label = '-1.0k Asteroids_-10k Points of Time Sampled'

fileLoc = root + extraBit + label + postfix

thing = tools.floatParseCSV(fileLoc)
print(len(thing))

thing = tools.transpose(thing)
print(len(thing))

timeList = thing[0] / (24 * 3600) / 365

title = (extraBit + label).replace('_', '\n')
title = title.replace('Acceleration', 'Magnitude')

magList = thing[1]**2 + thing[2]**2 + thing[3]**2
magList = np.sqrt(magList)

plt.plot(timeList, magList * 1e15, label = 'Absolute Magnitude')
# plt.plot(timeList, thing[1] * 1e15, label = 'x Acceleration')
# plt.plot(timeList, thing[2] * 1e15, label = 'y Acceleration')
# plt.plot(timeList, thing[3] * 1e15, label = 'z Acceleration')
plt.title(title)
plt.legend()
plt.ylim([0, 1500])
# plt.ylim([-1500, 1500])
# plt.ylim([-1e10, 1e10])
plt.xlabel('Time (years)')
plt.ylabel('Acceleration (fm/s^2)')
plt.show()
