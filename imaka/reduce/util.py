from astropy.table import Table
from astropy.io import fits
import os, errno, shutil
import pdb

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
