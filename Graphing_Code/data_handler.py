import my_module.util as tools
import numpy as np

_unitLoc = 'unit/'

def collect(location):
	data = tools.floatParseCSVfile(location)
	data = tools.transpose(data)
	return data

class classHandler:
	def __init__(self, fileName):
		self.fileName = fileName
		temp
