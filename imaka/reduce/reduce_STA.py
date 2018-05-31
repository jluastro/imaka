import numpy as np 
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.modeling import models
import ccdproc


def treat_overscan(files):
    """
    Performs bias/dark subtraction from overscan images on STA camera
    Assumes STA 8x2 'cell' structure, each with its own overscan region
    Works independent of image binning
    Input: an array of files from STA camera
    Output: writes new images with '_scan' suffix in name
    """
    for file in files:
        
        # Read in data and divide into 'cells'
        print('Working on file: ', file)
        dat, hdr = fits.getdata(file, header=True)
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
        scan_wid = int(cell_wid * (3/25))  # Width of overscan region (3/25 of cell)
        data_wid = int(cell_wid * (22/25)) # Width of data region (22/25 of cell)

        # For each cell, take the median of each overscan row and 
        # subtract that from the corresponding row in the data
        for ii in range(16):
            scan, img = np.split(cells[ii], [scan_wid], axis=1)
            scans.append(scan)
            imgs.append(img)
            med_col = np.median(scan, axis=1)
            med_cell = np.array([med_col,]*data_wid).transpose()
            clean_cell = img - med_cell 
            clean_cells.append(clean_cell)

        # Rebuild image with clean cells
        clean_one_row = np.hstack(clean_cells)
        clean_rows = np.hsplit(clean_one_row, 2)
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
        print("Working on file: ", file.split("/")[-1])
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
