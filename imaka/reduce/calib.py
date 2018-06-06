import os, sys
from . import util
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
import numpy as np
import pdb

module_dir = os.path.dirname(__file__)

def makedark(dark_files, output):
    """
    Make dark image for data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a darks/ subdirectory. 

    files: integer list of the files. Does not require padded zeros.
    output: output file name. Include the .fits extension.
    """
    # Catch if there are no files sent in.
    if len(dark_files) == 0:
        raise RuntimeError('No files passed into makedark')
        
    # Output directory is the same as the input directory.
    dark_dir = os.path.dirname(dark_files[0]) + '/'
    
    _out = output
    _outlis = output.replace('.fits', '.list')
    util.rmall([_out, _outlis])

    # Print out useful info.
    print('Making dark frame: ' + output)
    for ff in range(len(dark_files)):
        filedir, filename = os.path.split(dark_files[ff])
        print('\t' + filename)

    # Read in the dark images
    print('\nReading in files')
    darks = np.array([fits.getdata(dark_file) for dark_file in dark_files])

    # Sigma clip  (about the median)
    print('Sigma clipping and median combining (3 sigma lower, 2 sigma upper, 5 iterations) about the median.')
    dark_mean, dark_median, dark_std = sigma_clipped_stats(darks, sigma_lower=3, sigma_upper=2, iters=5, axis=0)
    dark = dark_median
    
    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(dark_files[0])
    fits.writeto(_out, dark.data, header=hdr, clobber=True)

    outlis_file = open(_outlis, 'w')
    for ff in range(len(dark_files)):
        outlis_file.write(dark_files[ff] + '\n')
    outlis_file.close()
    
    return

def makeflat(flat_files, dark_files, output_file, darks=True):
    """
    Make a dark-subtracted flat field image from a stack of flats
    and a stack of darks.

    This program should be run from the <epoch>/fli/ directory and
    it expects the "Pleiades" folder to be included in the file names.
    The output directory for the resulting flats are stored in

    Parameters
    ----------
    flat_files : list
        List of strings with full-path names of the individual flat-field exposures.
    dark_files : list
        List of strings with full-path names of the individual flat-field exposures.
        These should be matched 1:1 with the flats. 
    output_file : str
        Name (including directory) of the output flat field FITS file.
    """
    print('\nREDUCE_FLI: makeflat():')

    if darks==True:
        # Catch if there are no files sent in.
        if (len(flat_files) == 0) or (len(dark_files) == 0):
            raise RuntimeError('No files passed into makeflat')

    # Output File Names
    _out = output_file
    _outlis = output_file.replace('.fits', '.list')
    util.rmall([_out, _outlis])

    # Read in the dark images
    print('  Reading in dark files...')
    if darks==True:
        darks = np.array([fits.getdata(dark) for dark in dark_files])

    # Read in the flat images
    print('  Reading in flat files...')
    flats = np.array([fits.getdata(flat) for flat in flat_files])

    # Print out useful info.
    _lis = open(_outlis, 'w')
    _lis.write('#   Flat    Dark\n')

    if darks==True:
        for ff in range(len(flat_files)):
            filedir1, filename1 = os.path.split(dark_files[ff])
            filedir2, filename2 = os.path.split(flat_files[ff])
            print('  Flat: ' + filename2 + '  Dark: ' + filename1)
            _lis.write('{0:s}  {1:s}\n'.format(filename1, filename2))
            _lis.close()


    #Dark subtraction
    print('  Subtracting darks from flats...')
    if darks==True:
        flat_ds = flats - darks
    else:
        flat_ds = flats

    #Normalize each of the flats:
    print('  Normalizing each flat...')
    norm_flats = []
    for i in range(len(flat_ds)):
        mean, median, stdev = sigma_clipped_stats(flat_ds[i], sigma=3, iters=2, axis=None)
        norm = flat_ds[i] / median
        norm_flats.append(norm)

    # Sigma clip  (about the median)
    print('  Sigma clipping (3 sigma, 2 iterations) about the median...')
    mean, median, stdev = sigma_clipped_stats(norm_flats, sigma=3, iters=2, axis=0)

    # Median combine the masked images.
    print('  Median combining all images...')
    flat = median.data

    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(flat_files[0])
    fits.writeto(_out, flat, header=hdr, clobber=True)

    return


