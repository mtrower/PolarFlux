__authors__ = ["Zach Werginz", "Andres Munoz-Jaramillo"]
__email__ = ["zachary.werginz@snc.edu", "amunozj@gsu.edu"]

import numpy as np
import sunpy.map
from sunpy.sun import constants
from sunpy.sun import sun
import uncertainties.unumpy as unp
from uncertainties.umath import *
from uncertainties import ufloat
import astropy.units as u
import kpvt_class

class CRD:
    """Calculates various magnetogram coordinate information.
    Can calculate heliographic coordinate information,
    line of sight (LOS) corrections for the magnetic field,
    area elements for each pixel, and magnetic flux. This can
    be done for one pixel, or the whole data map. If the whole
    data map is given as a parameter, it will save the information
    as an instance attribute for the object.

    """

    RSUN_METERS = sun.constants.radius.si.value
    DSUN_METERS = sun.constants.au.si.value

    
    def __init__(self, filename):
        """Reads magnetogram as a sunpy.map object."""
        self.im_raw = sunpy.map.Map(filename)
        self.im_raw_u = np.array(np.abs(self.im_raw.data)*.05)

        if self.im_raw.detector == '512':

            # Define center of sun and location of detector.
            self.X0 = self.im_raw.meta['CRPIX1A']
            self.Y0 = self.im_raw.meta['CRPIX2A']
            self.B0 = self.im_raw.meta['B0']
            self.L0 = self.im_raw.meta['L0']
            self.xScale = self.im_raw.scale[0].value
            self.yScale = self.im_raw.scale[1].value
            self.rsun = self.im_raw.rsun_obs.value
            self.dsun = self.DSUN_METERS

        elif self.im_raw.detector == 'SPMG':

            # Define center of sun and location of detector.
            self.X0 = self.im_raw.meta['CRPIX1A']
            self.Y0 = self.im_raw.meta['CRPIX2A']
            self.B0 = self.im_raw.meta['B0']
            self.L0 = self.im_raw.meta['L0']
            self.xScale = self.im_raw.scale[0].value
            self.yScale = self.im_raw.scale[1].value
            self.rsun = self.im_raw.rsun_obs.value / self.im_raw.meta['SCALE']
            self.dsun = self.DSUN_METERS

        elif self.im_raw.detector == 'MDI':

            # Define center of sun and location of detector.
            self.X0 = self.im_raw.meta['X0']
            self.Y0 = self.im_raw.meta['Y0']
            self.B0 = self.im_raw.meta['B0']
            self.L0 = self.im_raw.meta['L0']
            self.xScale = 1.982
            self.yScale = 1.982
            self.rsun = self.im_raw.rsun_obs.value
            self.dsun = self.im_raw.dsun.value
            if self.im_raw.meta['p_angle'] == 180.0:
                self.im_raw.rotate(180)

        elif self.im_raw.detector == 'HMI':
            self.X0 = self.im_raw.meta['CRPIX1']
            self.Y0 = self.im_raw.meta['CRPIX2']
            self.B0 = self.im_raw.meta['CRLT_OBS']
            self.L0 = self.im_raw.meta['CRLN_OBS']
            self.xScale = self.im_raw.scale[0].value
            self.yScale = self.im_raw.scale[1].value
            self.rsun = self.im_raw.rsun_obs.value
            self.dsun = self.im_raw.dsun.value
            
    def __repr__(self):
        #TODO
        return None

    def meta(self):
        """Prints the fits header in a more readable fashion."""
        for key, value in self.im_raw.meta.items():
            print ("{}: {}".format(key, value))

    def heliographic(self, *args, array=True, corners=False):
        """Calculate hg coordinates from hpc and returns it.

        Can accept either a coordinate pair (x, y) or an entire 2D array.
        This pair corresponds the the pixel you want information on.

        The function will return a coordinate pair if given a coordinate
        pair or two arrays (latitude and longitude) if given a data array.
        Those arrays are indexed such that [0,0] is the top left pixel,
        although it refers to the bottom left portion of the picture.

        Use standard python indexing conventions for both the single
        coordinate and array calculations [row, column].

        Examples: 
        lath, lonh = kpvt.heliographic()
        aia.heliographic(320, 288, array=False)
        """

        # Check for single coordinate or ndarray object.
        if array:
            x, y = self._grid(corners)            
        else:
            # Coordinate conventions go [row, col] with
            # row zero being at the bottom (fits file)
            xScale = self.im_raw.scale[0].value
            yScale = self.im_raw.scale[1].value
            x = (args[1] - self.X0)*xScale
            y = (args[0] - self.Y0)*yScale

        # Calculations taken from sunpy.wcs.
        # First convert to heliocentric cartesian coordinates.
        rx, ry, rz = self._hpc_hcc(x, y)

        # Now convert to heliographic coordinates.
        lonh, lath = self._hcc_hg(rx, ry, rz)

        # Only add the instance attribute if it doesn't exist.
        if array and not hasattr(self, 'lonh') and not corners:
            self.lonh = lonh
            self.lath = lath
            return

        return lonh, lath

    def los_corr(self, *args, array=True):
        """Takes in coordinates and returns corrected magnetic field.

        Applies the dot product between the observers unit vector and
        the heliographic radial vector to get the true magnitude of 
        the magnetic field vector. See geometric projection for
        calulations.
        """

        print("Correcting line of sight magnetic field.")
        if array:
            try:
                lonh, lath = np.deg2rad(self.lonh), np.deg2rad(self.lath)
            except AttributeError:
                self.heliographic()
                lonh, lath = self.lonh*np.pi/180, self.lath*np.pi/180
        else:
            lonh, lath = np.deg2rad(self.heliographic(args[0], args[1]))

        B0 = ufloat(self.B0, np.abs(self.B0)*.05)*np.pi/180
        L0 = ufloat(self.L0, np.abs(self.L0)*.05)*np.pi/180

        Xobs = cos(B0)*cos(L0)
        Yobs = cos(B0)*sin(L0)
        Zobs = sin(B0)

        coslat = np.cos(lath)
        coslon = np.cos(lonh)
        sinlat = np.sin(lath)
        sinlon = np.sin(lonh)


        corr_factor = (coslat*coslon*Xobs.n + coslat*sinlon*Yobs.n + sinlat*Zobs.n)
        corr_factor_u = np.sqrt((coslat*coslon*Xobs.s)**2 + (coslat*sinlon*Yobs.s)**2 + (sinlat*Zobs.s)**2)
        if array:
            self.im_corr = self.im_raw.data/corr_factor
            self.im_corr_u = self.im_raw.
            bad_ind = np.where(self.rg > self.rsun*np.sin(75.0*np.pi/180))
            self.im_corr[bad_ind] = np.nan
            return
        else:
            return self.im_raw.data[args[0], args[1]]/corr_factor

    def eoa(self, *args, array=True):
        """Takes in coordinates and returns the area of pixels on sun.

        Each pixel is projected onto the sun, and therefore pixels close to
        the limbs have vastly greater areas. This function uses a closed form
        solution to a spherical area integral to calulate the area based on
        the heliographic coordinate unit vectors of each corner of the pixel.
        We use these to calculate a solid angle of a pyramid with its apex
        at the center of the sun.
        """

        print ("Calculating element of area.")
        # Assume coordinate is in center of pixel.
        # Information on pixel standard is in this article.
        # http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1061GFUL
        if array:
            lon, lat = np.deg2rad(self.heliographic(corners=True))
            # Calculating unit vectors of pixel corners for solid angle.
            r1 = self._spherical_to_cartesian(lon, lat, 0, 0)
            r2 = self._spherical_to_cartesian(lon, lat, 1, 0)
            r3 = self._spherical_to_cartesian(lon, lat, 1, 1)
            r4 = self._spherical_to_cartesian(lon, lat, 0, 1)

        else:
            x = args[0]
            y = args[1]
            lonUL, latUL = self.heliographic(x - .5, y - .5)
            lonLL, latLL = self.heliographic(x + .5, y - .5)
            lonLR, latLR = self.heliographic(x + .5, y + .5)
            lonUR, latUR = self.heliographic(x - .5, y + .5)

            # Calculating unit vectors of pixel corners for solid angle.
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

        # Calculate solid angle of pixel based on a pyrimid shaped polygon.
        # See http://planetmath.org/solidangleofrectangularpyramid

        cross1 = np.cross(r1, r2, axis=0)
        cross2 = np.cross(r3, r4, axis=0)
        numerator1 = self._dot(cross1, r3)
        numerator2 = self._dot(cross2, r1)
        solid_angle1 = 2*np.arctan2(numerator1,
                (self._dot(r1, r2) + self._dot(r2, r3) + self._dot(r3, r1) + 1))
        solid_angle2 = 2*np.arctan2(numerator2, 
                (self._dot(r3, r4) + self._dot(r4, r1) + self._dot(r3, r1) + 1))
        solid_angle = solid_angle1 + solid_angle2

        r = self.RSUN_METERS*100 # Convert to centimeters
        if array:
            self.area = np.abs((r**2)*solid_angle)
            ind = np.where(self.rg[1:len(self.rg)-1, 1:len(self.rg)-1] > self.rsun)
            self.area[ind] = np.nan
            return
        else:
            if self.rg > self.rsun:
                return np.nan
            else:    
                return np.abs((r**2)*solid_angle)

    def magnetic_flux(self, *args, array=True, raw_field=False):
        """Takes in coordinates and returns magnetic flux of pixel."""
        if array:
            try:
                area = self.area
            except AttributeError:
                self.eoa()
                area = self.area

            if raw_field:
                field = self.im_raw.data
                print ("Calculating raw magnetic flux.")
                self.mflux_raw = area*field
                return
            else:
                try:
                    field = self.im_corr
                except AttributeError:
                    self.los_corr()
                    field = self.im_corr
                finally:
                    print ("Calculating corrected magnetic flux.")
                    self.mflux_corr = area*field
                    return

        else:
            return self.eoa(*args)*self.los_corr(*args)
    
    def _grid(self, corners=False):
        """Create an xy grid of coordinates for heliographic array.

        Uses meshgrid. If corners is selected, this function will shift
        the array by half a pixel in both directions so that the corners
        of the normal array can be accessed easily.
        """
        # Retrieve integer dimensions and create arrays holding
        # x and y coordinates of each pixel
        xDim = np.int(np.floor(self.im_raw.dimensions[0].value))
        yDim = np.int(np.floor(self.im_raw.dimensions[1].value))
        
        if corners:
            xRow = (np.arange(0, xDim + 1) - self.X0 - 0.5)*self.xScale
            yRow = (np.arange(0, yDim + 1) - self.Y0 - 0.5)*self.yScale
        else:
            xRow = (np.arange(0, xDim) - self.X0)*self.xScale
            yRow = (np.arange(0, yDim) - self.Y0)*self.yScale
        
        self.xg, self.yg = np.meshgrid(xRow, yRow, indexing='xy')
        self.rg = np.sqrt(self.xg**2 + self.yg**2)

        return self.xg, self.yg

    def _hpc_hcc(self, x, y):
        """Converts hpc coordinates to hcc coordinates. 

        x -- x coordinate in arcseconds
        y -- y coordinate in arcseconds
        Calculations taken and shortened from sunpy.wcs.
        """
        x *= np.deg2rad(1)/3600.0
        y *= np.deg2rad(1)/3600.0

        q = self.dsun * np.cos(y) * np.cos(x)
        distance = q**2 - self.dsun**2 + self.RSUN_METERS**2
        
        distance = q - np.sqrt(distance)

        rx = distance * np.cos(y) * np.sin(x)
        ry = distance * np.sin(y)
        rz = np.sqrt(self.RSUN_METERS**2 - rx**2 - ry**2)

        return rx, ry, rz

    def _hcc_hg(self, x, y, z):
        """Converts hcc coordinates to Stonoyhurst heliographic. 

        x - x coordinate in meters
        y - y coordinate in meters 
        z - z coordinate in meters
        Calculations taken and shortened
        from sunpy.wcs.
        """
        cosb = np.cos(np.deg2rad(self.B0))
        sinb = np.sin(np.deg2rad(self.B0))

        hecr = np.sqrt(x**2 + y**2 + z**2)
        hgln = np.arctan2(x, z*cosb - y*sinb) \
                + np.deg2rad(self.L0)
        hglt = np.arcsin((y * cosb + z * sinb)/hecr)

        return hgln*180/np.pi, hglt*180/np.pi

    def _dot(self, a, b):
        """Vectorized version of the dot product of two arrays.

        Wanted a dot product function that performed operations
        element-wise.
        """
        return  np.array(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

    def _spherical_to_cartesian(self, lon, lat, i, j):
        """Takes latitude, longitude arrays and returns cartesian unit vector.

        Latitude and longitude must be in degrees and must have already been 
        shifted .5 pixels as calculated from corners keyword of heliographic
        function. i and j represent the shift and thus corner.
        i, j
        0, 0: top-left
        1, 0: bottom-left
        1, 1: bottom-right
        0, 1: top-right
        """
        coslat = np.cos(lat)
        coslon = np.cos(lon)
        sinlat = np.sin(lat)
        sinlon = np.sin(lon)
        l = len(lat)
        r = np.array([coslat[i:l - 1 + i, j:l - 1 + j]*coslon[i:l - 1 + i, j:l - 1 + j],
                    coslat[i:l - 1 + i, j:l - 1 + j]*sinlon[i:l - 1 + i, j:l - 1 + j],
                    sinlat[i:l - 1 + i, j:l - 1 + j]])
        return r