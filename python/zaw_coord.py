import numpy as np
import sunpy.map
import sunpy.wcs

import astropy.units as u

from sunpy.sun import constants
from sunpy.sun import sun

class CRD:
	
	def __init__(self, filename):
		self.im_raw = sunpy.map.Map(filename)
		

	def heliographic(self, x, y):
		
		# Must be in range
		if x < 0 or x > self.im_raw.meta['X0'] or y < 0 or y > self.im_raw.meta['Y0']:
			raise IndexError("Index out of range")

		#
		x = np.floor(x - self.im_raw.meta['X0'])
		y = np.floor(y - self.im_raw.meta['Y0'])
		return sunpy.wcs.convert_hpc_hg(x*self.im_raw.scale[0].value, y*self.im_raw.scale[1].value, b0_deg = self.im_raw.meta['B0'], l0_deg = self.im_raw.meta['l0'])

	def los_field(self, x, y):
		#Check for NaN, no calculation needed
		if (self.im_raw.data[x, y] == np.nan):
			return np.nan

		B0 = self.im_raw.meta['B0']
		L0 = self.im_raw.meta['L0']
		Xobs = np.cos(np.deg2rad(B0))*np.cos(np.deg2rad(L0))
		Yobs = np.cos(np.deg2rad(B0))*np.sin(np.deg2rad(L0))
		Zobs = np.sin(np.deg2rad(B0))

		Lonh = np.deg2rad(self.heliographic(x,y)[0])
		Lath = np.deg2rad(self.heliographic(x,y)[1])

		corr_factor = np.cos(Lath)*np.cos(Lonh)*Xobs + np.cos(Lath)*np.sin(Lonh)*Yobs + np.sin(Lath)*Zobs
		return self.im_raw.data[x, y]/corr_factor

