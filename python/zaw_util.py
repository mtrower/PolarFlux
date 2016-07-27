import os.path
import glob
import datetime as dt
from zaw_coord import CRD

__all__=['date2md', 'md2date', 'CRD_read', 'search_file']

data_root='.'

def date2md(date, instr):
	#Converts a standard date string into an instrument mission date.
	if instr == '512':
		DayOff = dt.date(1970, 1, 1)
	elif instr == 'spmg':
		DayOff = dt.date(1990, 1, 1)
	elif instr == 'mdi':
		DayOff = dt.date(1993, 1, 1)
	elif instr == 'hmi':
		DayOff = dt.date(2009, 1, 1)
	else:
		DayOff = dt.date(1970, 1, 1)
	
	return date.toordinal() - DayOff.toordinal()

def md2date(md, instr):
	#Converts a standard date string into an instrument mission date.
	if instr == '512':
		DayOff = dt.date(1970, 1, 1)
	elif instr == 'spmg':
		DayOff = dt.date(1990, 1, 1)
	elif instr == 'mdi':
		DayOff = dt.date(1993, 1, 1)
	elif instr == 'hmi':
		DayOff = dt.date(2009, 1, 1)
 	
	return dt.fromordinal(md + DayOff.toordinal())

def CRD_read(date, instr):
	try:
		filename = search_file(date, instr)[0]
	except IndexError:
		return -1

	mgnt = CRD(filename)
	if mgnt == -1:
		return -1
	mgnt.heliographic(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data, raw_field=True)
	mgnt.date = date
	mgnt.md = date2md(date, instr)

	return mgnt

def search_file(date, instr):
	if instr == '512':
		fn0 = data_root + '/KPVT'
		subdir = str(date.year - 1900) + str(date.month).zfill(2)

		return glob.glob(os.path.join(fn0, subdir, '*'+ str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '*.fits'))
	elif instr == 'spmg':
		fn0 = data_root + '/SPMG'
		subdir = datestr[2:4] + datestr[5:7]

		return glob.glob(os.path.join(fn0, subdir, '*'+ str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '*.fits'))
	# TODO figure out mdi_datestr nonsense
	elif instr == 'mdi':
		fn0 = data_root + '/MDI'
		subdir = ''

	elif instr == 'hmi':
		fn0 = data_root + '/HMI'

		return glob.glob(os.path.join(fn0, '*'+ str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '*.fits'))
	else:
		return -1
