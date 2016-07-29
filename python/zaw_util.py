import os.path
import glob
import datetime as dt
from astropy.io import fits
from zaw_coord import CRD

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
    # Set defaults
    subdir = ''
    fn0 = instr.upper()
    filename ='*%s*.fits' % date.strftime('%Y%m%d')

    # Set overrides
    if instr == '512':
        fn0 = 'KPVT'
        subdir = '%d%02d' % (date.year - 1900, date.month)
        filename = date.strftime('%Y%m%d') + '*.fits'

    elif instr == 'spmg':
        subdir = '%d%02d' % (date.year - 1900, date.month)

    elif instr == 'mdi':
        md = date2md(date, instr) + 1
        subdir = os.path.join(
                str(date.year)
                , 'fd_M_96m_01d.%06d' % md
        )
        filename ='fd_M_96m_01d.%d*.fits' % md
        
    elif instr == 'hmi':
        pass

    else:
        raise ValueError('Unrecognized instrument')

    # Execute
    searchspec = os.path.join(data_root, fn0, subdir, filename)
    files = glob.glob(searchspec)

    if not files:
        raise IOError('File not found')

    if instr == 'mdi':
        return mdi_file_choose(files)
    else:
        return files

def mdi_file_choose(f):
    best = None
    ival = 0
    mv = 100000
    for x in f:
        print (x)
        m = fits.open(x)
        try:
            intv = m[0].header['INTERVAL']
            if intv == '':
                intv = 0
            else:
                intv = int(intv)
            if intv >= ival:
                if int(m[0].header['MISSVALS']) < mv:
                    best = x
                    ival = m[0].header['INTERVAL']
                    mv = m[0].header['MISSVALS']
        except KeyError:
            continue
    if best == None:
        print(f[-1])
        return f[-1]
    return best
