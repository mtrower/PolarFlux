# This is a script that calculates flux information for magnetograms
# poles. The poles are based on certain criteria to make sure we
# don't factor in enormously projected areas of the sun.
#
#    Variables:
#				pc_px_ind:		Polar cap pixels
#				vpc_px_ind:		Valid polar cap pixels

import numpy as np
import pandas as pd
import datetime as dt
from zaw_coord import CRD

# This is the container class for the entire segment polar data.
class PFData:
	DayOff_512 = dt.date(1970, 1, 1)
	DayOff_SPMG = dt.date(1990, 1, 1)
	DayOff_MDI	= dt.date(1993, 1, 1)
	DayOff_HMI = dt.date(2009, 1, 1)
	pass

# Data dictionary initialization to be stored in a list.
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


def date(seg):
	if seg.instrument == '512':
		seg.DayOff = seg.DayOff_512
	elif seg.instrument == 'spmg':
		seg.DayOff = seg.DayOff_SPMG
	elif seg.instrument == 'mdi':
		seg.DayOff = seg.DayOff_MDI
	elif seg.instrument == 'hmi':
		seg.DayOff = seg.DayOff_HMI

	seg.mdi_i = seg.start_date.toordinal() - seg.DayOff.toordinal()
	seg.mdi_f = seg.end_date.toordinal() - seg.DayOff.toordinal()

	return

def read_file(date, instrument):
	# TODO : Use date to find filename

	try:
		mgnt = CRD(filename)
	except:
		return -1
	mgnt.magnetic_flux(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data, raw_field=True)
	mgnt.mdi_i = 
	mgnt.date = 

	return mgnt

def calc_pol(data, mgnt, pole):
	pole = pole.lower()
	swt = 0
	if pole == 'north':
		pc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > segment.deg_lim)

		vpc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > segment.deg_lim
						 and np.isfinite(mgnt.im_raw.data))

		pc_px_pos = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > segment.deg_lim 
						and mgnt.im_corr > 0.0)

		pc_px_neg = np.where(mgnt.rg < mgnt.rsun and mgnt.lath > segment.deg_lim 
						and mgnt.im_corr < 0.0)
	else:
		pc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -segment.deg_lim)

		vpc_px_ind = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -segment.deg_lim 
						and np.isfinite(mgnt.im_raw.data) )

		pc_px_pos = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -segment.deg_lim 
						and mgnt.im_corr > 0.0)

		pc_px_neg = np.where(mgnt.rg < mgnt.rsun and mgnt.lath < -segment.deg_lim 
						and mgnt.im_corr < 0.0)

	# Search for polarity mixture.
	if (np.size(pc_px_pos) == 0 or np.size(pc_px_neg) == 0):
		print ("%s hemisphere has no polarity mixture.", %pole)
		if pole == 'north':
			swt = 3
		else:
			swt = 3

	# Disregard magnetograms under the tolerance level for valid pixels.
	if (float(np.size(vpc_px_ind))/float(np.size(pc_px_ind)) < segment.inv_px_tol):
		print ("%s hemisphere has more than %s %% invalid pixels.", 
				%(pole, 1.0 - (segment.inv_px_tol)*100))
		if pole == 'north':
			swt = 2
		else:
			swt = 2

	# Return with no data if the two previous cases are true.
	if n_swt != 0 or s_swt != 0:
		data['mdi_i'] = mgnt.mdi_i
		data['date'] = mgnt.date
		return data
	# Otherwise calculate data from the poles.
	else:
		intf = np.nanmean(mgnt.im_raw[vpc_px_ind])
    	intfc = np.nanmean(mgnt.im_corr[vpc_px_ind])
    	unsflux = np.nansum(np.abs(mgnt.mflux_raw[vpc_px_indn]))
        unsfluxc = np.nansum(np.abs(mgnt.mflux_corr[vpc_pxindn]))
        sflux = np.nansum(mgnt.mflux_raw[vpc_px_ind])
        sfluxc = np.nansum(mgnt.mflux_corr[vpc_px_ind])
        posfluxc = np.nansum(mgnt.mflux_corr[pc_px_posn])
        negfluxc = np.nansum(mgnt.mflux_corr[pc_px_neg])
        visarea = np.nansum(mgnt.area[vpc_px_ind])
        max_pxflux = np.nanmax(np.abs( mgnt.mflux_corr[vpc_px_ind]))
        max_pxf = np.nanmax(mgnt.im_raw.data[vpc_px_ind])
        max_pxfc =  np.nanmax(mgnt.im_corr[vpc_px_ind])

        # Corrected field might be zero. Replace entries with NaNs.
        if unsfluxc eq 0.0:
        	unsfluxc = np.nan
        	sfluxc = np.nan
        	posfluxc = np.nan
        	negfluxc = np.nan
        	max_pxflux = np.nan
        	swt = 4

       	var = {'intf': intf, 'intfc': intfc, 'unsflux': unsflux,
       					'unsfluxc': unsfluxc, 'sflux': sflux,
       					'sfluxc': sfluxc, 'posfluxc': posfluxc,
       					'negfluxc': negfluxc, 'visarea': visarea,
       					'max_pxflux': max_pxflux, 'max_pxf': max_pxf,
       					'max_pxfc': max_pxfc}
       	# Assign data values per north/south pole.
   		for x in var.items():
   			if pole == 'north':
   				data[x[0] + '_n'] = x[1]
   			else:
   				data[x[0] + '_s'] = x[1]

   		return data


d1 = input('Enter starting date')
d2 = input('Enter end date')
instr = input('Enter the instrument')

# Initialize the segment data first.
segment = PFData()
segment.pf = []
segment.deg_lim = 65.0
segment.inv_px_tol = 0.85
segment.start_date = dt.date(int(d1[0:4]), int(d1[5:7]), int(d1[8:10]))
segment.end_date = dt.date(int(d2[0:4]), int(d2[5:7]), int(d2[8:10]))
segment.instrument = instr

# Proceed with function calls to read and calculate data.
date(segment)
md_c = segment.mdi_i

while md_c < seg.mdi_f:
	mgnt = read_file(md_c , segment.instrument)
	if mgnt != -1:
		segment.pf.append(calc_pol(data, mgnt, 'north'))
		segment.pf.append(calc_pol(data, mgnt, 'south'))
	md_current += 1


