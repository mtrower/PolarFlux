# This is a script that calculates flux information for magnetograms
# poles. The poles are based on certain criteria to make sure we
# don't factor in enormously projected areas of the sun.
#
#    Variables:
#				pc_px_ind:		Polar cap pixels
#				vpc_px_ind:		Valid polar cap pixels

import numpy as np
import pandas as pd
from zaw_coord import CRD

pf = []
deg_lim = 65.0
inv_px_tol = .85
mdi_i = 0
data = {'mdi_i': 0L, 'date': '', 'intf_n': np.nan, 'intfc_n': np.nan,
		'unsflux_n': np.nan, 'unsfluxc_n': np.nan, 'sflux_n': np.nan,
		'sfluxc_n': np.nan, 'posfluxc_n': np.nan, 'negfluxc_n': np.nan,
		'vnpc_px_n': np.nan, 'visarea_n': np.nan, 'max_pxflux_n': np.nan,
		'max_pxf_n': np.nan, 'max_pxfc_n': np.nan, 'n_swt': np.nan,
		'intf_s': np.nan, 'intfc_s': np.nan, 'unsflux_s': np.nan,
		'unsfluxc_s': np.nan, 'sflux_s': np.nan, 'sfluxc_s': np.nan,
		'posfluxc_s': np.nan, 'negfluxc_s': np.nan, 'vnpc_px_s': np.nan,
		'visarea_s': np.nan, 'max_pxflux_s': np.nan, 'max_pxf_s': np.nan,
		'max_pxfc_s': np.nan, 's_swt': np.nan}

def read_file(date):
	# TODO : Use date to find filename

	mgnt = CRD(filename)
	mgnt.magnetic_flux(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data, raw_field=True)
	mgnt.mdi_i = 
	mgnt.date = 

	return mgnt

def calc_pol(data, mgnt, pole):
	pole = pole.lower()
	swt = 0
	if pole == 'north':
		pc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > deg_lim)

		vpc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > deg_lim
						 and np.isfinite(mgnt.im_raw.data))

		pc_px_pos = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > deg_lim 
						and mgnt.im_corr > 0.0)

		pc_px_neg = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > deg_lim 
						and mgnt.im_corr < 0.0)
	else:
		pc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -deg_lim)

		vpc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -deg_lim 
						and np.isfinite(mgnt.im_raw.data) )

		pc_px_pos = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -deg_lim 
						and mgnt.im_corr > 0.0)

		pc_px_neg = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -deg_lim 
						and mgnt.im_corr < 0.0)

	if (np.size(pc_px_pos) == 0 or np.size(pc_px_neg) == 0):
		print ("%s hemisphere has no polarity mixture.", %pole)
		if pole == 'north':
			n_swt = 3
		else:
			s_swt = 3
	if (float(np.size(vpc_px_ind))/float(np.size(pc_px_ind)) < inv_px_tol):
		print ("%s hemisphere has more than %s %% invalid pixels.", %(pole, 1.0 - (inv_px_tol)*100))
		if pole == 'north':
			n_swt = 2
		else:
			s_swt = 2
	if n_swt != 0 or s_swt != 0:
		data['mdi_i'] = mgnt.mdi_i
		data['date'] = mgnt.date
		return data
