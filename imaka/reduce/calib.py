import os, sys
from . import util
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
import numpy as np
import multiprocessing as mp
import pdb
import itertools

module_dir = os.path.dirname(__file__)

try:
    str_type = types.StringTypes
    float_type = types.FloatType
    int_type = types.IntType
except:
    str_type = str
    float_type = float
    int_type = int

def makedark(dark_files, output, rescale=None):
    """
    Make dark image for data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a darks/ subdirectory. 

    Inputs
    ------
    files: list of ints
        integer list of the files. Does not require padded zeros.
    output: str
        output file name. Include the .fits extension.

    Optional Inputs
    --------------------
    rescale: None, float, or list of floats
        Useful for constructing skies from different integration times.
        Default is no rescaling.
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

    if rescale != None:
        if (isinstance(rescale, int_type) or isinstance(rescale, float_type)):
            print(f'Rescaling all frames a factor of {rescale:.1f}')
            darks *= rescale
        else:
            for ff in range(len(rescale)):
                print(f'Rescaling frame {ff} by a factor of {rescale[ff]:.1f}')
                darks[ff] *= rescale[ff]

    # Sigma clip  (about the median)
    print('Sigma clipping and median combining (3 sigma lower, 1 sigma upper, 5 iterations) about the median.')
    dark_mean, dark_median, dark_std = sigma_clipped_stats(darks, sigma_lower=3, sigma_upper=1, maxiters=5, axis=0)
    dark = dark_mean
    
    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(dark_files[0])
    fits.writeto(_out, dark.data, header=hdr, overwrite=True)

    outlis_file = open(_outlis, 'w')
    for ff in range(len(dark_files)):
        outlis_file.write(dark_files[ff] + '\n')
    outlis_file.close()
    
    return

def makeflat(flat_files, dark_files, output_file, darks=True, fourfilter=False):
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
        dark_imgs = np.array([fits.getdata(dark) for dark in dark_files])

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
        flat_ds = flats - dark_imgs
    else:
        flat_ds = flats

    # Prepare for parallel compute.
    # Use N_cpu - 2 so we leave 2 for normal
    # operations. 
    cpu_count = _cpu_num()

    # Start the pool
    pool = mp.Pool(cpu_count)

    # Normalize each of the flats: in parallel
    print(f'  Normalizing each flat in parallel with {cpu_count} cores.')
    if fourfilter:
        norm_flats = pool.map(normalize_flat_filt, [flat_ds[i] for i in range(len(flat_ds))])
    else:
        norm_flats = pool.map(normalize_flat, [flat_ds[i] for i in range(len(flat_ds))])

    # Sigma clip  (about the median)
    print('  Sigma clipping (3 sigma, 2 iterations) about the median...')
    mean, median, stdev = sigma_clipped_stats(norm_flats, sigma=3, maxiters=2, axis=0)

    # Median combine the masked images.
    print('  Median combining all images...')
    flat = median.data

    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(flat_files[0])
    fits.writeto(_out, flat, header=hdr, overwrite=True)

    return


def make_mask(flat_file, out_filename, mask_min=None, mask_max=None,
              left_slice=None, right_slice=None, top_slice=None, bottom_slice=None):
    """
    Make a mask of all the regions we shouldn't be looking for stars. 
    Be sure to run this before find_stars. 
    """
    flat, hdr = fits.getdata(flat_file, header=True)

    mask = mask_pix(flat, mask_min, mask_max, 
                    left_slice=left_slice, right_slice=right_slice,
                    top_slice=top_slice, bottom_slice=bottom_slice)

    fits.writeto(out_filename, mask.astype(np.uint8), header=hdr, overwrite=True)

    return

def mask_pix(flat_img, mask_min, mask_max, left_slice=0, right_slice=0, top_slice=0, bottom_slice=0):
    """
    Masks pixels in flat image based on input minimum and maximum values.
    Also masks edges as specified by 'slice' values. Note that indexing/orientation 
    is based in python, not fits (0,0 is top left, not bottom left, and starts at 0)
    
    Returns boolean mask array
    """

    # Read in file
    y, x = np.shape(flat_img)

    # Mask pixels below a defined threshold 'mask_min'
    copy1 = np.copy(flat_img)
    copy1[np.where((copy1>mask_min) & (copy1<mask_max))] = 0
    mask1 = np.ma.make_mask(copy1,copy=True, shrink=True,dtype=np.bool)

    # Mask edge pixels as defined by 'slice' variables 
    copy2 = np.copy(flat_img)
    copy2[top_slice:y-bottom_slice, left_slice:x-right_slice] = 0
    mask2 = np.ma.make_mask(copy2, copy=True, shrink=True, dtype=np.bool)

    # Combine masks
    combo_mask = np.ma.mask_or(mask1, mask2)
    
    return combo_mask


def normalize_flat(flat_img):
    mean, median, stdev = sigma_clipped_stats(flat_img, sigma=3, maxiters=2, axis=None)
    norm = flat_img / median
    return norm

def normalize_flat_filt(flat_img):
    # Normalize by quadrant instead of whole image
    ## TODO: smaller box for normalization
    
    med_f = np.zeros_like(flat_img)
    
    # Define quadrants with indicies
    img_half = int(flat_img.shape[0]/ 2)
    q_edge = 1000 #0 if use all of quad in median
    
    # Creating a list of masks for Quadrants 
    ind_1 = np.zeros_like(flat_img)
    ind_1[0:img_half, 0:img_half] = 1
    
    inds = [(0, img_half), (img_half, img_half*2)]
    
    for i in itertools.product(inds, inds):
        #smaller mask for stats
        ind_mask = np.zeros_like(flat_img)
        ind_mask[i[0][0]+q_edge:i[0][1]-q_edge, i[1][0]+q_edge:i[1][1]-q_edge] = 1
        mask = ~np.array(ind_mask, dtype=bool)
        mean, median, stdev = sigma_clipped_stats(flat_img, sigma=3, maxiters=2, axis=None, mask=mask)
        
        #larger mask for median
        ind = np.zeros_like(flat_img)
        ind[i[0][0]:i[0][1], i[1][0]:i[1][1]] = 1
        med_f += ind*median
        
    norm = flat_img / med_f
    return norm

def combine_filter_flat(flat_long, flat_short, output_file, filt_order, flip_180=True):
    """
    Combine flat field image for differing filters. Overwrite one quad of short flat with "I" quad from long flat

    Parameters
    ----------
    flat_long : str
        String of full path for make_flat output, longer exposure for I band
    flat_short : str
        String of full path for make_flat output, shorter exposure
    output_file : str
        Name (including directory) of the output flat field FITS file.
    filt_order : str
        BVIR order of filter, used to determine I quadrant
    """
    print('\nREDUCE_FLI: combine_filter_flat():')
    
    # Output File Names
    _out = output_file
    _outlis = output_file.replace('.fits', '.list')
    util.rmall([_out, _outlis])
    
    # Read in the flat images
    print('  Reading in flat files...')
    flat_long_data = np.array(fits.getdata(flat_long))
    flat_short_data = np.array(fits.getdata(flat_short))
    
    if flat_long_data.shape != flat_short_data.shape:
        print("ERROR: flats not the same size")
    
    # Define quadrants with indicies
    img_half = int(flat_short_data.shape[0]/ 2)
    
    # Iterate over all filters to add to main flat
    for i, filt in enumerate(filt_order):
        if filt == "I":
            row = (i//2 + 1) % 2
            col = (i+row + 1) % 2 
            if flip_180:
                row = (row+1)%2
                col = (col+1)%2
            r = [img_half*row, img_half*(row+1)]
            c = [img_half*col, img_half*(col+1)]
            flat_short_data[r[0]:r[1], c[0]:c[1]] = flat_long_data[r[0]:r[1], c[0]:c[1]]
        else:
            continue
    
    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(flat_long)
    fits.writeto(_out, flat_short_data, header=hdr, overwrite=True)
    
    return

def _cpu_num():
    # How many cores to use for parallel functions
    # returns a number
    # Use N_cpu - 2 so we leave 2 for normal
    # operations. 
    cpu_count = mp.cpu_count()
    if (cpu_count > 2):
        cpu_count -= 2
        cpu_count = cpu_count if cpu_count <= 8 else 8
        
    # reurn 1 #(DEBUG)
    return cpu_count