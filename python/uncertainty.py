import numpy as np

__all__ = []

class uval:
	def __init__(self, value, uncertainty):
		self.v = value
		self.u = uncertainty

def add(*args):
	result = 0
	for term in args:
		result += term.uncertainty**2

	return np.sqrt(result)

def mult(*args):
