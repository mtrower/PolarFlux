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
	
		return sunpy.wcs.convert_hpc_hg(x*self.im_raw.scale[0].value, y*self.im_raw.scale[1].value, b0_deg = self.im_raw.meta['B0'], l0_deg = self.im_raw.meta['l0'])