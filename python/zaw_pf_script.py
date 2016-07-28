# This is a script that calculates flux information for magnetograms
# poles. The poles are based on certain criteria to make sure we
# don't factor in enormously projected areas of the sun.
#
#    Variables:
#               pc_px_ind:      Polar cap pixels
#               vpc_px_ind:     Valid polar cap pixels

import numpy as np
import pandas as pd
import datetime as dt
import zaw_util
from collections import OrderedDict
import os.path
import getopt
import sys

d1 = None       # Start date
d2 = None       # End date
instr = None    # Instrument

# This is the container class for the entire segmented polar data.
class PFData:
    pass

# Data dictionary initialization to be stored in a list.
data0 = OrderedDict([('md', 0), ('date', ''), ('intf_n', np.nan),
        ('intfc_n', np.nan), ('unsflux_n', np.nan), ('unsfluxc_n', np.nan),
        ('sflux_n', np.nan), ('sfluxc_n', np.nan), ('posfluxc_n', np.nan),
        ('negfluxc_n', np.nan), ('nvp_px_n', np.nan), ('visarea_n', np.nan), 
        ('max_pxflux_n', np.nan), ('max_pxf_n', np.nan), ('max_pxfc_n', np.nan),
        ('swt_n', np.nan), ('intf_s', np.nan), ('intfc_s', np.nan),
        ('unsflux_s', np.nan), ('unsfluxc_s', np.nan), ('sflux_s', np.nan),
        ('sfluxc_s', np.nan), ('posfluxc_s', np.nan), ('negfluxc_s', np.nan),
        ('nvp_px_s', np.nan), ('visarea_s', np.nan), ('max_pxflux_s', np.nan), 
        ('max_pxf_s', np.nan), ('max_pxfc_s', np.nan), ('swt_s', np.nan)])

def calc_pol(pf_data, mgnt, pole):
    pole = pole.lower()

    p_px, vp_px, posp_px, negp_px = indices(mgnt, pole, segment.deg_lim)
    
    swt = validate(p_px, vp_px, posp_px, negp_px)
    
    # Return with no data if the two previous cases are true.
    if swt != 0:
        pf_data['md'] = mgnt.md
        pf_data['date'] = mgnt.date
        if pole == 'north':
            pf_data['swt_n'] = swt
        else:
            pf_data['swt_s'] = swt

        return pf_data
    # Otherwise calculate data from the poles.
    else:
        intf = np.nanmean(mgnt.im_raw.data[vp_px])
        intfc = np.nanmean(mgnt.im_corr[vp_px])
        unsflux = np.nansum(np.abs(mgnt.mflux_raw[vp_px]))
        unsfluxc = np.nansum(np.abs(mgnt.mflux_corr[vp_px]))
        sflux = np.nansum(mgnt.mflux_raw[vp_px])
        sfluxc = np.nansum(mgnt.mflux_corr[vp_px])
        posfluxc = np.nansum(mgnt.mflux_corr[posp_px])
        negfluxc = np.nansum(mgnt.mflux_corr[negp_px])
        visarea = np.nansum(mgnt.area[vp_px])
        max_pxflux = np.nanmax(np.abs( mgnt.mflux_corr[vp_px]))
        max_pxf = np.nanmax(mgnt.im_raw.data[vp_px])
        max_pxfc =  np.nanmax(mgnt.im_corr[vp_px])

        # Corrected field might be zero. Replace entries with NaNs.
        if unsfluxc == 0.0:
            unsfluxc = np.nan
            sfluxc = np.nan
            posfluxc = np.nan
            negfluxc = np.nan
            max_pxflux = np.nan
            swt = 4

        var = {'intf': intf, 'intfc': intfc, 'unsflux': unsflux,
                        'unsfluxc': unsfluxc, 'sflux': sflux,
                        'sfluxc': sfluxc, 'posfluxc': posfluxc,
                        'negfluxc': negfluxc, 'nvp_px': np.size(vp_px),
                        'visarea': visarea, 'max_pxflux': max_pxflux,
                        'max_pxf': max_pxf, 'max_pxfc': max_pxfc, 'swt': swt}
        # Assign data values per north/south pole.
        for x in var.items():
            if pole == 'north':
                pf_data[x[0] + '_n'] = x[1]
            else:
                pf_data[x[0] + '_s'] = x[1]
        pf_data['md'] = mgnt.md
        pf_data['date'] = mgnt.date     
        return pf_data

def indices(m, pole, dlim):
    if pole == 'north':
        pc = np.where(np.logical_and(np.less(m.rg, m.rsun),
                np.greater(m.lath, dlim)))

        vpc = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.greater(m.lath, dlim)),
                np.isfinite(m.im_raw.data)))

        pc_pos = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.greater(m.lath, dlim)),
                np.greater(m.im_corr, 0.0)))

        pc_neg = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.greater(m.lath, dlim)), 
                np.less(m.im_corr, 0.0)))
    else:
        pc = np.where(np.logical_and(np.less(m.rg, m.rsun),
                np.less(m.lath, -dlim)))

        vpc = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.less(m.lath, -dlim)),
                np.isfinite(m.im_raw.data)))

        pc_pos = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.less(m.lath, -dlim)),
                np.greater(m.im_corr, 0.0)))

        pc_neg = np.where(np.logical_and(
                np.logical_and(np.less(m.rg, m.rsun), np.less(m.lath, -dlim)), 
                np.less(m.im_corr, 0.0)))
    return pc, vpc, pc_pos, pc_neg

def validate(p, vp, pos, neg):
    # Search for polarity mixture.
    if (np.size(pos) == 0 or np.size(neg) == 0):
        print ("{} hemisphere has no polarity mixture.".format(pole))
        return 3

    # Disregard magnetograms under the tolerance level for valid pixels.
    if (float(np.size(vp))/float(np.size(p)) < segment.inv_px_tol):
        print ("{} hemisphere has more than {} %% invalid pixels.". 
                format(pole, 1.0 - (segment.inv_px_tol)*100))
        return 2

    return 0

def usage():
    print('Usage: zaw_pf_script.py [-d data-root] [-s start-date] [-e end-date] [-i instrument]')

def parse_args():
    global d1, d2, instr

    try:
        opts, args = getopt.getopt(
                sys.argv[1:],
                "d:s:e:i:", 
                ["data-root=", "date-start=", "date-end=", "instrument="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-d", "--data-root"):
            zaw_util.data_root = arg
        elif opt in ("-s", "--date-start"):
            d1 = arg
        elif opt in ("-e", "--date-end"):
            d2 = arg
        elif opt in ("-i", "--instrument"):
            instr = arg
        else:
            assert False, "unhandled option"

parse_args()

if    d1 == None:    d1 = input('Enter starting date: ')
if    d2 == None:    d2 = input('Enter end date: ')
if instr == None: instr = input('Enter the instrument: ')

# Initialize the segment data first.
segment = PFData()
segment.pf = []
segment.deg_lim = 65.0
segment.inv_px_tol = 0.85
segment.start_date = dt.date(int(d1[0:4]), int(d1[5:7]), int(d1[8:10]))
segment.end_date = dt.date(int(d2[0:4]), int(d2[5:7]), int(d2[8:10]))
segment.instrument = instr
segment.md_i = zaw_util.date2md(segment.start_date, segment.instrument)
segment.md_f = zaw_util.date2md(segment.end_date, segment.instrument)

# Proceed with function calls to read and calculate data.
date_c = segment.start_date

while date_c <= segment.end_date:
    print(date_c)
    mgnt = zaw_util.CRD_read(date_c, segment.instrument)
    if mgnt != -1:
        print("Calculating polar parameters")
        data = data0.copy()
        calc_pol(data, mgnt, 'north')
        calc_pol(data, mgnt, 'south')
        segment.pf.append(data)
    date_c = date_c + dt.timedelta(1)

obj = pd.DataFrame(segment.pf, columns=segment.pf[0].keys())
out_fn = 'PF_{}_{}.csv'.format(segment.start_date.isoformat(), segment.end_date.isoformat())
obj.to_csv(out_fn)