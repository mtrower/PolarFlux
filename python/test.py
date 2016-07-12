from zaw_coord import CRD
import sunpy.wcs
import numpy as np
import timeit



x = CRD('fd_M_96m_01d.1222.0005.fits')

#Heliographic testing
# start = timeit.default_timer()
#lath, lonh = x.heliographic(x.im_raw.data)
#print ( "Array Latitude = %s Longitude = %s " %(lath[512,512], lonh[511,511]) )
#lonh, lath = x.Heliographic(0, 0)
#print ( "Coordinate Longitude = %s Latitude = %s " %(lonh, lath) )


# print( sunpy.wcs.convert_hpc_hg(0*x.im_raw.scale[0].value, 0*x.im_raw.scale[1].value, b0_deg = x.im_raw.meta['B0'], l0_deg = x.im_raw.meta['L0']) )
# stop = timeit.default_timer()

#LOS testing
# corr = x.los_corr(x.im_raw.data)

# print( "Raw field = %s " %x.im_raw.data[750,750])
# print( "Corrected field = %s " %corr[750, 750])

# print(corr)

#Element of Area testing
# areapix = x.area(238,238)
# print("%e" %areapix.value)
# areamap = x.area(x.im_raw.data)
# print("%e" %areamap[749,749].value)

# Magnetic flux testing
flux = x.magnetic_flux(x.im_raw.data)
print(flux[900,900])
print(flux.nansum())

#print ("Time = ", stop - start)