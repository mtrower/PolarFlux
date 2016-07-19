__authors__ = ["Zach Werginz", "Andres Munoz-Jaramillo"]
__email__ = ["zachary.werginz@snc.edu", "amunozj@gsu.edu"]

import numpy as np
import sunpy.map
from sunpy.sun import constants
from sunpy.sun import sun
import astropy.units as u

#np.set_printoptions(threshold=np.nan)

RSUN_METERS = sun.constants.radius.si.value
DSUN_METERS = sun.constants.au.si.value


class CRD:
	"""Calculates various magnetogram coordinate information.



	"""
	
	def __init__(self, filename):
		"""Reads magnetogram as a sunpy.map object."""
		self.im_raw = sunpy.map.Map(filename)
		try:
			self.B0 = self.im_raw.meta['B0']
		except KeyError:
			self.B0 = self.im_raw.meta['OBS_B0']
		try:
			self.L0 = self.im_raw.meta['L0']
		except KeyError:
			self.L0 = self.im_raw.meta['OBS_L0']
		try:
			self.X0 = self.im_raw.meta['X0']
		except KeyError:
			self.X0 = self.im_raw.meta['IMG_X0']
		try:
			self.Y0 = self.im_raw.meta['Y0']
		except KeyError:
			self.Y0 = self.im_raw.meta['IMG_Y0']
			
		self.rsun = self.im_raw.rsun_obs.value
		

	def __repr__(self):
		#TODO
		return None

	def heliographic(self, *args):
		"""Calculate hg coordinates from hpc and returns it.

		Can accept either a coordinate pair (x, y) or an entire 2D array.

		The cartesian coordinates are measured from the center pixel
		of the map.
		The function will return a coordinate pair if given a coordinate
		pair or two arrays (latitude and longitude) if given a data array.
		Those arrays are indexed such that [0,0] is the top left pixel.

		Examples: 
		lath, lonh = kpvt.heliographic(kpvt.im_raw.data)
		aia.heliographic(320, 288)
		"""

		xScl = self.im_raw.scale[0].value
		yScl = self.im_raw.scale[1].value
		
		# Check for single coordinate or ndarray object.
		if isinstance(args[0], np.ndarray):
			#Retrieve integer dimensions and create arrays holding
			#x and y coordinates of each pixel
			xdim = np.int(np.floor(self.im_raw.dimensions[0].value))
			ydim = np.int(np.floor(self.im_raw.dimensions[1].value))
			try: #if len(args) == 1:
				xrow = (np.arange(0, xdim) - self.X0 + args[1])*xScl
				yrow = (np.arange(0, ydim) - self.Y0 + args[2])*yScl
				self.xg, self.yg = np.meshgrid(xrow, yrow, indexing='xy')
				self.rg = np.sqrt(self.xg**2 + self.yg**2)
				x = self.xg
				y = -self.yg
			except IndexError: #Don't include the shift
				xrow = (np.arange(0, xdim) - self.X0)*xScl
				yrow = (np.arange(0, ydim) - self.Y0)*yScl
				self.xg, self.yg = np.meshgrid(xrow, yrow, indexing='xy')
				self.rg = np.sqrt(self.xg**2 + self.yg**2)
				x = self.xg
				y = -self.yg
		else:
			x = (args[1] - self.X0)*xScl
			y = (self.Y0 - args[0])*yScl
		
		# First convert to heliocentric cartesian coordinates.
		cosx = np.cos(x*np.deg2rad(1)/3600.0)
		sinx = np.sin(x*np.deg2rad(1)/3600.0)
		cosy = np.cos(y*np.deg2rad(1)/3600.0)
		siny = np.sin(y*np.deg2rad(1)/3600.0)

		q = DSUN_METERS * cosy * cosx
		distance = q**2 - DSUN_METERS**2 + RSUN_METERS**2
		distance = q - np.sqrt(distance)

		rx = distance * cosy * sinx
		ry = distance * siny
		rz = np.sqrt(RSUN_METERS**2 - rx**2 - ry**2)

		# Now convert to heliographic coordinates.
		cosb = np.cos(np.deg2rad(self.B0))
		sinb = np.sin(np.deg2rad(self.B0))

		hecr = np.sqrt(rx**2 + ry**2 + rz**2)
		hgln = np.arctan2(rx, rz*cosb - ry*sinb) \
		       + np.deg2rad(self.L0)
		hglt = np.arcsin((ry * cosb + rz * sinb)/hecr)
		# Only add the instance attribute if it doesn't exist.
		if isinstance(args[0], np.ndarray) and not hasattr(self, 'lonh'):
			self.lonh = np.rad2deg(hgln)
			self.lath = np.rad2deg(hglt)

		return np.rad2deg(hgln), np.rad2deg(hglt)

	def los_corr(self, *args):
		"""Takes in coordinates and returns corrected magnetic field.

		Applies the dot product between the observers unit vector and
		the heliographic radial vector to get the true magnitude of 
		the magnetic field vector. See geometric projection.
		"""

		if isinstance(args[0], np.ndarray):
			lonh, lath = np.deg2rad(self.heliographic(args[0]))
		else:
			lonh, lath = np.deg2rad(self.heliographic(args[0], args[1]))

		Xobs = np.cos(np.deg2rad(B0))*np.cos(np.deg2rad(self.L0))
		Yobs = np.cos(np.deg2rad(B0))*np.sin(np.deg2rad(self.L0))
		Zobs = np.sin(np.deg2rad(B0))

		corr_factor = (np.cos(lath)*np.cos(lonh)*Xobs
					   + np.cos(lath)*np.sin(lonh)*Yobs
					   + np.sin(lath)*Zobs)
		if isinstance(args[0], np.ndarray):
			self.im_corr = self.im_raw.data/corr_factor
			return self.im_corr
		else:
			return self.im_raw.data[args[0], args[1]]/corr_factor

	def eoa(self, *args):
		""" Takes in coordinates and returns the area of pixels on sun.

		Each pixel is projected onto the sun, and therefore pixels close to
		the limbs have vastly greater areas. This function uses a closed form
		solution to a spherical area integral to calulate the area based on
		the heliographic coordinates of the lower-left(LL) part of each pixel
		and upper-right(UR) part of each pixel.
		"""

		#Assume coordinate is in center of pixel.
		#Information on pixel standard is in this article.
		#http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1061GFUL
		if isinstance(args[0], np.ndarray):
			lonUL, latUL = self.heliographic(args[0], -.5, -.5)
			lonLL, latLL = self.heliographic(args[0], .5, -.5)
			lonLR, latLR = self.heliographic(args[0], .5, .5)
			lonUR, latUR = self.heliographic(args[0], -.5, .5)
		else:
			x = args[0]
			y = args[1]
			lonUL, latUL = self.heliographic(x - .5, y - .5)
			lonLL, latLL = self.heliographic(x + .5, y - .5)
			lonLR, latLR = self.heliographic(x + .5, y + .5)
			lonUR, latUR = self.heliographic(x - .5, y + .5)

		r1 = np.array([np.cos(np.deg2rad(latUL))*np.cos(np.deg2rad(lonUL)),
						np.cos(np.deg2rad(latUL))*np.sin(np.deg2rad(lonUL)),
						np.sin(np.deg2rad(latUL))])

		r2 = np.array([np.cos(np.deg2rad(latLL))*np.cos(np.deg2rad(lonLL)),
						np.cos(np.deg2rad(latLL))*np.sin(np.deg2rad(lonLL)),
						np.sin(np.deg2rad(latLL))])

		r3 = np.array([np.cos(np.deg2rad(latLR))*np.cos(np.deg2rad(lonLR)),
						np.cos(np.deg2rad(latLR))*np.sin(np.deg2rad(lonLR)),
						np.sin(np.deg2rad(latLR))])

		r4 = np.array([np.cos(np.deg2rad(latUR))*np.cos(np.deg2rad(lonUR)),
						np.cos(np.deg2rad(latUR))*np.sin(np.deg2rad(lonUR)),
						np.sin(np.deg2rad(latUR))])

		
		crosspn = np.cross(r1, r2, axisa=0, axisb=0, axisc=0)
		print("Dim cross product is", crosspn.shape)

		numerator1 = dot(crosspn, r3)
		numerator2 = dot(np.cross(r3, r4, axisa=0, axisb=0, axisc=0), r1)
		#result1 = numerator1/(dot(r1, r2) + dot(r2, r3) + dot(r3, r1) + 1)
		#result2 = numerator2/(dot(r2, r3) + dot(r3, r4) + dot(r4, r2) + 1)
		
		solid_angle1 = 2*np.arctan2(numerator1, (dot(r1, r2) + dot(r2, r3) + dot(r3, r1) + 1))
		solid_angle2 = 2*np.arctan2(numerator2, (dot(r3, r4) + dot(r4, r1) + dot(r3, r1) + 1))
		solid_angle = solid_angle1 + solid_angle2
		r = 6.957e10 * u.cm
		
		#dPhi is the change in azimuthal angle.
		dPhi = lonUR - lonLL
		dPhi = np.deg2rad(dPhi)
		#Theta is used for the change in polar angle.
		theta2 = np.deg2rad(90-latLL)
		theta1 = np.deg2rad(90-latUR)

		#area element = r^2sin(theta)*dtheta*dphi
		if isinstance(args[0], np.ndarray):
			self.area = np.abs((r**2)*solid_angle)
			#self.area = np.abs((r**2)*dPhi*(np.cos(theta1) - np.cos(theta2)))
			#self.area = (r**2)*2*np.arctan(np.sin(.5*(latUR + latLL))/np.sin(.5*(latUR - latLL))*np.tan((lonUR - lonLL)/5))
			#r**2*np.abs(np.sin(np.deg2rad(latUR)) - np.sin(np.deg2rad(latLL)))*np.abs(np.deg2rad(lonUR) - np.deg2rad(lonLL))
			ind = np.where(self.rg > self.rsun)
			self.area[ind] = np.nan
			return self.area
		else:
			return np.abs((r**2)*solid_angle)

	def magnetic_flux(self, *args):
		""" Takes in coordinates and returns magnetic flux of pixel."""
		area = self.area(*args)
		field = self.los_corr(*args)
		if isinstance(args[0], np.ndarray):
			self.mgnt_flux = area*field
		return area*field

def dot(a, b):
	return  np.array(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])