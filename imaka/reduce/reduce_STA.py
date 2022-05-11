import numpy as np 
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.modeling import models
import os
import math
import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import pdb
#from astroscrappy import detect_cosmics
import glob
import photutils
from photutils import psf
from photutils import morphology as morph
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
#from flystar import match
#from flystar import align
#from flystar import transforms
from imaka.reduce import calib
from imaka.reduce import util
#import ccdproc
from scipy.ndimage import interpolation
from scipy.ndimage import median_filter
import scipy.ndimage
from astropy.table import Table
from skimage.measure import block_reduce
from datetime import datetime
import pytz
from imaka.reduce import reduce_fli
from astropy.table import Column
import multiprocessing as mp
from itertools import repeat

def treat_overscan(files, remake=False):
    """
    Performs bias/dark subtraction from overscan images on STA camera
    Assumes STA 8x2 'cell' structure, each with its own overscan region
    Works independent of image binning
    Input: an array of files from STA camera
    Output: writes new images with '_scan' suffix in name
    """
    # Parallael process the images.

    # Use N_cpu - 2 so we leave 2 for normal
    # operations. 
    cpu_count = mp.cpu_count()
    if (cpu_count > 2):
        cpu_count -= 2

    # Start the pool
    pool = mp.Pool(cpu_count)

    # Run
    print(f'Treating overscan in parallel with {cpu_count} cores.')
    pool.starmap(treat_overscan_single, zip(files, repeat(remake)))

    pool.close()
    pool.join()
    
    return


# Define a function to treat a single frame.
# Need this for parallelization.
def treat_overscan_single(filename, remake):
    pid = mp.current_process().pid
    
    # Write new image
    new_filename = filename[:-5] + '_scan.fits'
    if filename.endswith('.gz'):
        new_filename = filename[:-8] + '_scan.fits'

    # Check to see if it exists. If not, make it.
    if os.path.exists(new_filename) and not remake:
        return
        
    # Read in data
    print(f'p{pid} - treat_overscan_single on file: {filename}')
    dat, hdr = fits.getdata(filename, header=True)

    # Determine binning from image shape
    binfac = util.get_bin_factor(dat, hdr)

    # Split the data vertically for the two different read-outs. Take the 
    # bottom 8 and gang them on the end of the top 8 for easier access.
    bottom_ylo = int(0)
    bottom_yhi = int(5280 / binfac)
    top_ylo = int(5320 / binfac)
    top_yhi = int(10600 / binfac)
    row_bottom = dat[bottom_ylo:bottom_yhi, :]
    row_top = dat[top_ylo:top_yhi, :]
    one_row = np.hstack([row_bottom, row_top])
    cells = np.hsplit(one_row, 16)

    # Define the 'image' and 'overscan' region 
    # Done with fractions, so independent of binning
    imgs = []
    scans = []
    clean_cells = []

    # Check the date and adjust parameters accordingly.
    year = int(hdr['DATE-OBS'].split('-')[0])
    if year <= 2020:
        orig_wid = np.shape(dat)[1]        # Width of original image
        cell_wid = int(orig_wid / 8)       # Width of one cell (1/8 of image)
        scan_wid = int(120 / binfac)       # Width of overscan region (3/25 of cell)
        data_wid = int(1440 / binfac)      # Width of data region (22/25 of cell)
        pix_shift = 5
    else:
        orig_wid = np.shape(dat)[1]        # Width of original image
        cell_wid = int(orig_wid / 8)       # Width of one cell (1/8 of image)
        scan_wid = int(240 / binfac)       # Width of overscan region (2/13 of cell)
        data_wid = int(1320 / binfac)      # Width of data region (11/13 of cell)
        pix_shift = 0
        
    print(f'  p{pid} - orig_wid = {orig_wid}, cell_wid = {cell_wid}, ' +
          f'data_wid = {data_wid}, pix_shift = {pix_shift}')

    # For each cell, take the median of each overscan row and 
    # subtract that from the corresponding row in the data
    for ii in range(16):
        ## Take each data column located between ni*dw+ow and (ni+1)*dw and put 
        ## them next to each other in the new image
        img_lo = ii*(data_wid + scan_wid) + pix_shift
        img_hi = img_lo + data_wid
        scn_lo = img_hi
        if scn_lo < 0:
            scn_lo = 0
        scn_hi = scn_lo + scan_wid

        img = one_row[:, img_lo:img_hi]
        scn = one_row[:, scn_lo:scn_hi]

        med_col = np.median(scn, axis=1)
        med_cell = np.array([med_col,]*img.shape[1]).transpose()

        clean_cell = img - med_cell
        clean_cells.append(clean_cell)

    # Rebuild image with clean cells
    clean_one_row = np.hstack(clean_cells)
    clean_rows = np.hsplit(clean_one_row, 2)
    clean_image = np.vstack(clean_rows)

    fits.writeto(new_filename, clean_image, hdr, overwrite=True)

    return


def treat_overscan_2020_edit(files):
    """
    Performs bias/dark subtraction from overscan images on STA camera
    Assumes STA 8x2 'cell' structure, each with its own overscan region
    Works independent of image binning
    Input: an array of files from STA camera
    Output: writes new images with '_scan' suffix in name
    """
    for file in files:

        # Read in data
        print('Working on file: ', file)
        dat, hdr = fits.getdata(file, header=True)

        # Determine binning from image shape
        binfac = util.get_bin_factor(dat, hdr)
        #pdb.set_trace()

        # Split the data vertically for the two different read-outs. Take the 
        # bottom 8 and gang them on the end of the top 8 for easier access.
        rows = np.vsplit(dat, 2)
        one_row = np.hstack([rows[0], rows[1]])
        cells = np.hsplit(one_row, 16)

        # Define the 'image' and 'overscan' region 
        # Done with fractions, so independent of binning
        imgs = []
        scans = []
        clean_cells = []

        orig_wid = np.shape(dat)[1]        # Width of original image
        cell_wid = orig_wid / 8            # Width of one cell (1/8 of image)
        scan_wid = int(120 / binfac)       # Width of overscan region (3/25 of cell)
        data_wid = int(1440 / binfac)      # Width of data region (22/25 of cell)
        pix_shift = 5
        print('orig_wid = ', orig_wid)
        print('cell_wid = ', cell_wid)
        print('scan_wid = ', scan_wid)
        print('data_wid = ', data_wid)

        # For each cell, take the median of each overscan row and 
        # subtract that from the corresponding row in the data
        for ii in range(16):
            ## Take each data column located between ni*dw+ow and (ni+1)*dw and put 
            ## them next to each other in the new image
            #img_lo = ii*data_wid + scan_wid - pix_shift
            img_lo = ii*(data_wid + scan_wid) + pix_shift
            #img_hi = (ii + 1)*data_wid - pix_shift
            img_hi = img_lo + data_wid
            #scn_lo = ii*data_wid - pix_shift
            scn_lo = img_hi
            if scn_lo < 0:
                scn_lo = 0
            scn_hi = scn_lo + scan_wid

            #img = one_row[img_lo:img_hi, :]
            #scn = one_row[scn_lo:scn_hi, :]
            img = one_row[:, img_lo:img_hi]
            scn = one_row[:, scn_lo:scn_hi]

            med_col = np.median(scn, axis=1)
            #med_cell = np.array([med_col,]*img.shape[0]).transpose()
            med_cell = np.array([med_col,]*img.shape[1]).transpose()

            clean_cell = img - med_cell
            clean_cells.append(clean_cell)

            #pdb.set_trace()
            # scan, img = np.split(cells[ii], [scan_wid - pix_shift], axis=1)
            # scans.append(scan)
            # imgs.append(img)
            # med_col = np.median(scan, axis=1)
            # med_cell = np.array([med_col,]*data_wid).transpose()
            # clean_cell = img - med_cell 
            # clean_cells.append(clean_cell)

        # Rebuild image with clean cells
        clean_one_row = np.hstack(clean_cells)
        pdb.set_trace()
        clean_rows = np.hsplit(clean_one_row, 2)
        clean_image = np.vstack(clean_rows)

        # Write new image
        new_filename = file[:-5] + '_scan.fits'
        fits.writeto(new_filename, clean_image, hdr, overwrite=True)
    
    return

def treat_overscan_working(files):
    """
    Performs bias/dark subtraction from overscan images on STA camera.
    Assumes STA 8x2 'cell' structure, each with its own overscan region.
    Works independent of image binning and checks image creation date
    before applying overscan treatment.
    Input: an array of files from STA camera
    Output: writes new images with '_scan' suffix in name
    """
    
    sta_date = datetime.strptime('2020-01-21','%Y-%m-%d')  #date of first known sky imgs taken with
                                                           #the 5280x5760 STA cam settings 
    for file in files:
        
        # Read in data and divide into 'cells'
        print('Working on file: ', file)
        dat, hdr = fits.getdata(file, header=True)
        date = datetime.strptime(hdr['DATE-OBS'],'%Y-%m-%d')
        
        rows = np.vsplit(dat, 2)
        one_row = np.hstack([rows[0], rows[1]])
        cells = np.hsplit(one_row, 16)

        if date >= sta_date:
            # Define the 'image' and 'overscan' region 
            # Done with fractions, so independent of binning
            # print("Using post-2018 STA overscan protocol...")
            imgs = []
            scans = []
            clean_cells = []
            corr_term = 24

            orig_wid = np.shape(dat)[1]             # Width of original image
            cell_wid = int(orig_wid / 8)            # Width of one cell (1/8 of image)
            scan_wid = int(round(cell_wid * (3/25))) - corr_term # Width of overscan region (3/25 of cell)
            data_wid = int(round(cell_wid * (22/25))) + corr_term # Width of data region (22/25 of cell)
            

            # For each cell, split the overscan and data regions 
            # into separate variables
            for ii in range(16):
                cell = cells[ii]
                scan = cell[:,:scan_wid]
                scans.append(scan)
                img = cell[:,scan_wid:data_wid]
                clean_cells.append(img)
        else:
            print("Using pre-2020 STA overscan protocol...")
            # Define the 'image' and 'overscan' region 
            # Done with fractions, so independent of binning
            imgs = []
            scans = []
            clean_cells = []

            orig_wid = np.shape(dat)[1]        # Width of original image
            cell_wid = orig_wid / 8            # Width of one cell (1/8 of image)
            scan_wid = int(round(cell_wid * (3/25)))  # Width of overscan region (3/25 of cell)
            data_wid = int(round(cell_wid * (22/25))) # Width of data region (22/25 of cell)

            # For each cell, take the median of each overscan row and 
            # subtract that from the corresponding row in the data
            for ii in range(16):
                scan, img = np.split(cells[ii], [scan_wid], axis=1)
                scans.append(scan)
                imgs.append(img)
                med_col = np.median(scan, axis=1)
                med_cell = np.array([med_col,]*(data_wid)).transpose()
                clean_cell = img - med_cell 
                clean_cells.append(clean_cell)
            
        # Rebuild image with clean cells
        clean_one_row = np.hstack(clean_cells)
        clean_rows = np.hsplit(clean_one_row,2)
        clean_image = np.vstack(clean_rows)
        
        # Write new image
        new_filename = file[:-5] + '_scan.fits'
        fits.writeto(new_filename, clean_image, hdr, overwrite=True)
    
    return


def clean_image(image_file):
    """
    Cleans a single image.  Cleaning includes:
        -Gain application
        -Cosmic ray removal
    
    Returns: astropy.nddata.ccddata.CCDData object
    
    Header needs: 'GAINI' (optional, but good); 'EXPTIME', 'BINFAC'
    """

    # Read in image and header
    img, hdr = fits.getdata(image_file, header=True)
    if 'GAINI' in hdr:
        gain = hdr['GAINI'] * u.electron/u.adu
    else:
        print("Gain not found in header: assuming 1 e/adu")
        gain = 1 *u.electron/u.adu

    readnoise = 5 * u.electron #From STA1600LN data sheet

    # Define data and uncertainty
    data = CCDData(img, unit=u.adu)
    data_with_deviation = ccdproc.create_deviation(data, gain=gain, readnoise=readnoise)
    data_with_deviation.header['exposure'] = float(hdr['EXPTIME'])

    # Apply gain to data
    gain_corrected = ccdproc.gain_correct(data_with_deviation, gain)

    # Clean cosmic rays
    cr_cleaned = ccdproc.cosmicray_lacosmic(gain_corrected, sigclip=5)

    return cr_cleaned


def create_bias(bias_files):
    """
    Cleans and median combines an input list of bias files with sigma clipping 3-σ about median
    
    Input: a list of bias image fits files 

    Returns: Saves master bias frame as 'master_bias.fits' in same directory as bias frames
    """
    
    # Read in bias images
    combiner = []
    for file in bias_files:
        print("create_bias on file: ", file.split("/")[-1])
        ccd = clean_image(file)
        combiner.append(ccd)

    # Sigma clip 3-σ about median
    print("Sigma clipping")
    combiner.sigma_clipping(func=np.ma.median)

    # Combine frames
    print("Combining frames")
    combined_median = combiner.median_combine()

    # Write to fits file with header info of last file
    print("Writing file")
    bias_dir = file.split('bias_')[0]
    combined_median.write(bias_dir+'master_bias.fits')

    return



def add_focus(stats_files):
    """
    Retrieves focus value from image headers for
    each frame in a stats_file, and creates a 
    new column with these values.  For STA camera.
    """
    for file in stats_files:
        tab = Table.read(file)
        N_frames = len(tab)
        focci = []
        for ii in range(N_frames):
            frame_file = tab['Image'][ii]
            hdr = fits.getheader(frame_file)
            focus = hdr['FOCUS']
            focci.append(focus)
        col_focus = Column(name='Focus', data=focci)
        tab.add_column(col_focus)
        tab.write(file, overwrite=True)

    return 


def four_filt_split(starlists, filt_order):
    
    """
    For STA camera with four filters.
    Splits starlist into four starlists,
    each corresponding to one quadrant/filter
    Inputs: starlists: array of STA camera starlists;
            filt_order: a four character string
            giving filter labels from in order:
            NW, NE, SE, SW, eg: 'BVIR'
    Output: Four new starlists for each input
            starlist given with ending:
            R_BVIR_stars.txt' for R filter and
            filter order 'BVIR'
    """

    filt_1 = filt_order[0]
    filt_2 = filt_order[1]
    filt_3 = filt_order[2]
    filt_4 = filt_order[3]

    N = len(starlists)
    for ii in range(N):
        print("Splitting file", ii+1, "of", N)
        print("\t", starlists[ii])

        # Read in image and starlists
        stars = Table.read(starlists[ii], format='ascii')
        x = stars['xcentroid']
        y = stars['ycentroid']
        stars.add_column(Column(np.repeat('4', len(stars)), name='FILTER'))

        # Define quadrants with indicies
        img_half = 5280 / 2 # FOR BIN2!!!!
        img_half = 5280 #BIN1

        ind_1 = np.where((x<img_half) & (y>img_half))
        ind_2 = np.where((x>img_half) & (y>img_half))
        ind_3 = np.where((x>img_half) & (y<img_half))
        ind_4 = np.where((x<img_half) & (y<img_half))
        
        stars['filter'][ind_1] = filt_1
        stars['filter'][ind_2] = filt_2
        stars['filter'][ind_3] = filt_3
        stars['filter'][ind_4] = filt_4
    
        # Write new files
        list_root = starlists[ii].split('stars.txt')[0]
        #stars[ind_1].write(list_root + filt_1 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        #stars[ind_2].write(list_root + filt_2 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        #stars[ind_3].write(list_root + filt_3 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        #stars[ind_4].write(list_root + filt_4 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        ## Retyring with more aggressive ascii formatting
        formats = {'xcentroid': '%8.3f', 'ycentroid': '%8.3f', 'sharpness': '%.2f',
                   'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
                   'flux': '%10.6f', 'mag': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
                   'theta': '%6.3f'}
        stars[ind_1].write(list_root + filt_1 + '_' + filt_order + '_stars.txt', format='ascii.fixed_width',
                          delimiter=None, bookend=False, formats=formats, overwrite=True)
        stars[ind_2].write(list_root + filt_2 + '_' + filt_order + '_stars.txt', format='ascii.fixed_width',
                          delimiter=None, bookend=False, formats=formats, overwrite=True)
        stars[ind_3].write(list_root + filt_3 + '_' + filt_order + '_stars.txt', format='ascii.fixed_width',
                          delimiter=None, bookend=False,  formats=formats, overwrite=True)
        stars[ind_4].write(list_root + filt_4 + '_' + filt_order + '_stars.txt', format='ascii.fixed_width',
                          delimiter=None, bookend=False,  formats=formats, overwrite=True)
        
        #from the original star function 
        #sources.write(img_file.replace('.fits', '_stars.txt'), format='ascii.fixed_width',
        #                  delimiter=None, bookend=False, formats=formats, overwrite=True)
    
    return
