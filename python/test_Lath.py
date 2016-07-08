from zaw_coord import CRD
import numpy as np

x = CRD('fd_M_96m_01d.1222.0005.fits')

x.Latarr = []

for i in np.nditer(x.im_raw.data):
	x.Latarr.append(x.heliographic(i%1024, i/1024))

x.Latarr = np.array(x.Latarr)
print(x)