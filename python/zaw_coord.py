import numpy as np
import sunpy.map

import astropy.units as u

from sunpy.sun import constants
from sunpy.sun import sun

rsun_meters = sun.constants.radius.si.value
dsun_meters = sun.constants.au.si.value

class CRD:
	
	def __init__(self, filename):
		self.im_raw = sunpy.map.Map(filename)
		

	def heliographic(self, *args):
		xScl = self.im_raw.scale[0].value
		yScl = self.im_raw.scale[1].value

		if isinstance(args[0], np.ndarray):
			#Retrieve integer dimensions and create arrays holding x and y coordinates of each pixel
			try: #if len(args) == 1:
				xdim = np.int(np.floor(self.im_raw.dimensions[0].value))
				ydim = np.int(np.floor(self.im_raw.dimensions[1].value))
				row = np.arange(0, xdim)
				x = ((np.tile(row, (xdim, 1)) - self.im_raw.meta['X0'] + args[1])*xScl).T
				y = (np.tile(row, (ydim, 1)) - self.im_raw.meta['Y0'] + args[1])*yScl
			except IndexError: #include the shift
				xdim = np.int(np.floor(self.im_raw.dimensions[0].value))
				ydim = np.int(np.floor(self.im_raw.dimensions[1].value))
				row = np.arange(0, xdim)
				x = ((np.tile(row, (xdim, 1)) - self.im_raw.meta['X0'])*xScl).T
				y = (np.tile(row, (ydim, 1)) - self.im_raw.meta['Y0'])*yScl
		else:
			x = (args[0] - self.im_raw.meta['X0']%np.floor(self.im_raw.meta['X0']))*xScl
			y = (args[1] - self.im_raw.meta['Y0']%np.floor(self.im_raw.meta['Y0']))*yScl
		
		#First convert to heliocentric cartesian coordinates
		cosx = np.cos(x * np.deg2rad(1) / 3600.0)
		sinx = np.sin(x * np.deg2rad(1) / 3600.0)
		cosy = np.cos(y * np.deg2rad(1) / 3600.0)
		siny = np.sin(y * np.deg2rad(1) / 3600.0)

		q = dsun_meters * cosy * cosx
		distance = q ** 2 - dsun_meters ** 2 + rsun_meters ** 2
		distance = q - np.sqrt(distance)

		rx = distance * cosy * sinx
		ry = distance * siny
		rz = np.sqrt(rsun_meters**2 - rx**2 - ry**2)

		#Now convert to heliographic coordinates
		cosb = np.cos(np.deg2rad(self.im_raw.meta['B0']))
		sinb = np.sin(np.deg2rad(self.im_raw.meta['B0']))

		hecr = np.sqrt(rx**2 + ry**2 + rz**2)
		hgln = np.arctan2(rx, rz * cosb - ry * sinb) + np.deg2rad(self.im_raw.meta['L0'])
		hglt = np.arcsin((ry * cosb + rz * sinb) / hecr)

		return np.rad2deg(hgln), np.rad2deg(hglt)

	def los_corr(self, *args):
		if isinstance(args[0], np.ndarray):
			lonh, lath = np.deg2rad(self.heliographic(args[0]))
		else:
			lonh, lath = np.deg2rad(self.heliographic(args[0], args[1]))

		#Check for NaN, no calculation needed
		#if (self.im_raw.data[x, y] == np.nan):
		#	return np.nan

		B0 = self.im_raw.meta['B0']
		L0 = self.im_raw.meta['L0']
		Xobs = np.cos(np.deg2rad(B0))*np.cos(np.deg2rad(L0))
		Yobs = np.cos(np.deg2rad(B0))*np.sin(np.deg2rad(L0))
		Zobs = np.sin(np.deg2rad(B0))

		corr_factor = np.cos(lath)*np.cos(lonh)*Xobs + np.cos(lath)*np.sin(lonh)*Yobs + np.sin(lath)*Zobs
		if len(args) == 2:
			return self.im_raw.data[args[0], args[1]]/corr_factor
		else:
			return self.im_raw.data/corr_factor

	def area(self, *args):

		#Assume coordinate is in center of pixel.
		#http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1061GFUL
		if isinstance(args[0], np.ndarray):
			lonLL, latLL = self.heliographic(args[0], -0.5)
			lonUR, latUR = self.heliographic(args[0], 0.5)
		else:
			x = args[0]
			y = args[1]
			lonLL, latLL = self.heliographic(x - .5, y - .5)
			lonUR, latUR = self.heliographic(x + .5, y + .5)

		r = 6.957e10 * u.cm
		dPhi = -(np.deg2rad(lonLL) - np.deg2rad(lonUR))
		theta2 = np.deg2rad(90+latLL)
		theta1 = np.deg2rad(90+latUR)


		return (r**2)*dPhi*(np.cos(theta1) - np.cos(theta2))

	def magnetic_flux(self, *args):
		area = self.area(*args)
		field = self.los_corr(*args)
		return area*field