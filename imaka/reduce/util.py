import numpy as np
from astropy.table import Table
from astropy.io import fits
import os, errno, shutil
import pdb
import pytz
from datetime import datetime

# Load up directory aliases
module_dir = os.path.dirname(__file__)
    
def rmall(files):
    """Remove list of files without confirmation."""
    for file in files:
        if os.access(file, os.F_OK): os.remove(file)
            
    return

def mkdir(dir):
    """Make directory if it doesn't already exist."""
    try: 
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

    return

def getcwd():
    """
    IRAF doesn't like long file names. This reduces them.
    """
    curdir = os.getcwd()
    
    curdir += '/'

    return curdir


def cp_change_prefix(arg1,arg2):
    """
    Takes files beginning with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory beginning with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files=[filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        pre = files[ff][0:ln]
        if pre == arg1:
            suf = files[ff][len(arg1):]
            newFile = arg2 + suf
            shutil.copy(files[ff], newFile)

    return

def cp_change_suffix(arg1,arg2):
    """
    Takes files ending with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory ending with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files=[filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        suf = files[ff][len(files[ff])-len(arg1):]
        if suf == arg1:
            pre = files[ff][0:len(files[ff])-len(arg1)]
            newFile = pre + arg2 
            shutil.copy(files[ff], newFile)

    return

def update_header_coords(fileList):
    """
    Updates coordinates in the header for XREF, YREF
    and XSTREHL, and YSTREHL.
 
    fileList : list of files to update
    """

    _files = Table.read(fileList, format='ascii', header_start=None)
    cols = _files.columns.keys()
    files = _files[cols[0]]
    files = [files[ff].split('.')[0] for ff in range(len(files))]
    

    for ff in range(len(files)):
        # Open .coo file and read 16C's coordinates
        coo = Table.read(files[ff]+'.coo', format='ascii', header_start=None)
        coo_cols = coo.columns.keys()
        xref = coo[coo_cols[0]]
        yref = coo[coo_cols[1]]

        # Open .coord file and read strehl source's coordinates
        coord = Table.read(files[ff]+'.coord', format='ascii', header_start=None)
        coord_cols = coord.columns.keys()
        xstr = coord[coord_cols[0]]
        ystr = coord[coord_cols[1]]
 
        # Open image and write reference star x,y to fits header
        fits_f = fits.open(files[ff]+'.fits')

        fits_f[0].header.update('XREF', "%.3f" %xref,
                              'Cross Corr Reference Src x')
        fits_f[0].header.update('YREF', "%.3f" %yref,
                              'Cross Corr Reference Src y')
        fits_f[0].header.update('XSTREHL', "%.3f" %xstr,
                              'Strehl Reference Src x')
        fits_f[0].header.update('YSTREHL', "%.3f" %ystr,
                              'Strehl Reference Src y')

        # Output fits file
        _out = 'new_hdr/' + files[ff] + '.fits'
        fits_f[0].writeto(_out, output_verify='silentfix')

    return

def get_plate_scale(img, hdr):
    """
    Return plate scale in arcsec / pixel. Should work 
    for all of our cameras. 
    """
    scale = 1.0
    
    if 'BINFAC' in hdr:
        scale_orig = 0.04 # " / pixel

        scale = scale_orig * hdr['BINFAC']

    # STA Camera
    if 'CCDBIN1' in hdr:
        scale = 0.063 * hdr['CCDBIN1']
    
    return scale

def get_times(hdr):
    camera = get_camera(hdr)

    if camera == 'STA':
        # Make dates and times in UT and HST
        if 'SYS-TIME' in hdr:
            time_hst = hdr['SYS-TIME']
        else:
            time_hst = hdr['SYSTIME']
    
        if 'SYS-DATE' in hdr:
            date_hst = hdr['SYS-DATE']
        else:
            date_hst = hdr['SYSDATE']
            
        date_utc = hdr['DATE-OBS']
    
        time = hdr['UT']
        h, m, s = time.split(':')
        new_time = h+":"+m+":"+s[:2]
        time_utc = new_time

    else:
        # Make dates and times in UT.
        # The time comes in from the header with "hh:mm:ss PM" in HST.
        # The date comes in from the header in HST.
        time_tmp = hdr['TIMEOBS']
        date_tmp = hdr['DATEOBS']
        hst_tz = pytz.timezone('US/Hawaii')

        if hdr['SHUTTER'] == True:
            dt_hst = datetime.strptime(date_tmp + ' ' + time_tmp, '%m/%d/%Y %H:%M:%S')
            dt_hst = hst_tz.localize(dt_hst)
        else:
            dt_hst = datetime.strptime(date_tmp + ' ' + time_tmp, '%d/%m/%Y %I:%M:%S %p')
            dt_hst = hst_tz.localize(dt_hst)
            noon = datetime.time(12, 0, 0) #assuming you'll never be taking images at local noon...
            del_day = datetime.timedelta(days=1)
            if dt_hst.time() < noon:
                dt_hst += del_day
                    
        dt_utc = dt_hst.astimezone(pytz.utc)


        time_hst = str(dt_hst.time())
        date_hst = str(dt_hst.date())
        time_utc = str(dt_utc.time())
        date_utc = str(dt_utc.date())

    return date_hst, time_hst, date_utc, time_utc

def get_camera(hdr):
    if 'TIMEOBS' in hdr:
        # This is an FLI only keyword.
        camera = 'FLI'
    else:
        camera = 'STA'

    return camera

def get_four_filter_quadrant(starlist_name):
    # Get filter position in case of four filter data
    name_strings = starlist_name.split("_")
    filt_name    = name_strings[-3]
    filt_order   = name_strings[-2]

    if filt_name == filt_order[0]:
        quad = 'NW'
    elif filt_name == filt_order[1]:
        quad = 'NE'
    elif filt_name == filt_order[2]:
        quad = 'SE'
    elif filt_name == filt_order[3]:
        quad = 'SW'

    return quad

def get_four_filter_name(starlist_name):
    # Get filter position in case of four filter data
    name_strings = starlist_name.split("_")
    filt_name    = name_strings[-3]
    filt_order   = name_strings[-2]

    return filt_name, filt_order
    
def get_filter(hdr):
    if 'FILT' in hdr:
        return hdr['FILT']

    if 'FILTER' in hdr:
        return hdr['FILTER']
    else:
        return 'I'
    
    return

def get_wavelength(filter_str):
    wave = 0
    if filter_str == 'B':
        wave = 445
    elif filter_str == 'V':
        wave = 551
    elif filter_str == 'R':
        wave = 658
    elif filter_str == 'I':
        wave = 806
    return wave

def get_quad(filt, filt_order):
    quad = ''
    if filt == filt_order[0]:
        quad = 'NW'
    elif filt == filt_order[1]:
        quad = 'NE'
    elif filt == filt_order[2]:
        quad = 'SE'
    elif filt == filt_order[3]:
        quad = 'SW'
    return quad
    
def get_bin_factor(image, header):

    if 'BINFAC' in header:
        return header['BINFAC']

    if 'CCDBIN1' in header:
        return header['CCDBIN1']
    
    imgshape = image.shape

    ## Determine binning from image size
    ## For 1x1 binning
    if (imgshape[0] > 10500) and (imgshape[1] > 11500):
        binfac = 1

    ## For 2x2 binning
    elif (imgshape[0] > 5200) and (imgshape[0] < 5500):
        binfac = 2
    
    ## For 3x3 binning
    elif (imgshape[0] > 3500) and (imgshape[0] < 3700):
        binfac = 3

    ## For 4x4 binning
    elif (imgshape[0] > 2600) and (imgshape[0] < 2700):
        binfac = 4

    return binfac

