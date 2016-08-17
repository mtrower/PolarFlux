from zaw_coord import CRD
import kpvt_class
#import sunpy.wcs
import numpy as np
from uncertainty import Measurement as M
import timeit
import datetime as dt
start = timeit.default_timer()

mdi = CRD('MDI\\fd_M_96m_01d.1222.0005.fits')
#kpvt = CRD('512c_eo000_C1_19771001_2048.fits')
#spmg = CRD('spmg_eo100_C1_19920421_1700.fits')
#hmi = CRD('HMI\\hmi.M_720s.20100504_214800_TAI.1.magnetogram.fits')

#MDI
#mdi.heliographic()
#print(mdi.lonh[996,506])
#mdi.eoa()
#lonh, lath = np.deg2rad(mdi.lonh), np.deg2rad(mdi.lath)
#mdi.los_corr()
#print(mdi.im_corr_u[996,506])
#mdi.magnetic_flux()
#print(mdi.mflux_corr[300,300])
#print(M.nansum(mdi.area))
#print (M.nansum(mdi.mflux_corr))
#print(np.nansum(mdi.mflux_corr[np.where(np.isfinite(mdi.mflux_corr))]))
#print (np.nansum(mdi.area))

#KPVT
#kpvt.heliographic()
#kpvt.los_corr() 
#kpvt.eoa()
#kpvt.magnetic_flux()
#print(np.nanmax(kpvt.area))
#print(kpvt.area.shape)
#print(np.nansum(kpvt.area))
#print(np.nanmax(kpvt.mflux_corr))

#SPMG
#spmg.heliographic()
#spmg.los_corr()
#spmg.eoa()
# spmg.magnetic_flux()
#HMI
#hmi.heliographic()
#hmi.los_corr()
#hmi.eoa()
#hmi.magnetic_flux()
#print(np.nansum(hmi.area))

#PF SCRIPT TESTING
import zaw_pf_script as p
import zaw_util
mgnt = CRD('MDI\\fd_M_96m_01d.1222.0005.fits')
mgnt.heliographic()
mgnt.magnetic_flux()
mgnt.magnetic_flux(raw_field=True)
p.d1 = '1996-05-07'
p.d2 = '1996-05-10'
x = p.PFData()
p.init(x)
mgnt.date = x.meta['start_date']
mgnt.md = zaw_util.date2md(mgnt.date, x.meta['instrument'])
p_px, vp_px, posp_px, negp_px = p.indices(mgnt, 'north', x.meta['deg_lim'])
data = p.data0.copy()
p.calc_pol(data, mgnt, 'north', x)


stop = timeit.default_timer()
print ("Time = ", stop - start)