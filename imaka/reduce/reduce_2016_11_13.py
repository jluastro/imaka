import pylab as plt
from PIL import Image
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
import photutils
from photutils import psf
from photutils import morphology as morph
from photutils import DAOStarFinder
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clipped_stats
import poppy
import pdb
from imaka.reduce import reduce_fli

from imaka.reduce import calib

def make_darks():
    root_dir = '/Users/jlu/GoogleDrive/Instruments/imaka/imaka/Commissioning/'
    root_dir += '2016-11 Observing/20161115/Focal_plane_images/obs_11152016/darks/'

    darks = glob.glob(root_dir + 'dark0*.fits')
    
    calib.makedark(darks[0:6], 'dark.fits')

    return


def reduce_pleiadies_east():
    data_dir = '/Users/jlu/Google Drive/Instruments/imaka/imaka (1)/Commissioning/2016-11 Observing/20161115/'
    data_dir += 'Focal_plane_images/obs_11152016/'

    dark_dir = data_dir + 'darks/'
    plei_dir = data_dir + 'Pleides_E/'

    
    dark_file = 'off_0S_0E.fits'


    # img_files = ['off_0S_0E_1.fits', 'off_0S_0E_2.fits', 'off_10S_0E_1.fits', 'off_10S_0E_2.fits',
    #              'off_10S_30E_1.fits', 'off_10S_30E_2.fits']
    img_files = ['off_0S_0E_1.fits']

    from scipy.ndimage import median_filter

    for ii in range(len(img_files)):
        print('Working on image: ', img_files[ii])
        # img_ds = reduce_fli.subtract_dark(img_files[ii], dark_file)
        img_ds, hdr_ds = fits.getdata(img_files[ii], header=True)

        # Fix hot and dead pixels
        print('\t Fix bad pixels')
        bad_pixels, img_cr = find_outlier_pixels(img_ds, tolerance=5, worry_about_edges=False)
 
        # Save the image
        fits.writeto(img_files[ii].replace('.fits', '_cr.fits'), img_cr, hdr_ds, clobber=True)
        fits.writeto(img_files[ii].replace('.fits', '_mask.fits'), bad_pixels, hdr_ds, clobber=True)

        # Bin by a factor of 10.
        print('\t Bin by a factor of 10')
        from skimage.measure import block_reduce
        img_bin = block_reduce(img_cr, (10, 10), func=np.sum)
        fits.writeto(img_files[ii].replace('.fits', '_bin.fits'), img_bin, hdr_ds, clobber=True)

        # Subtract background.
        print('\t Subtract background')
        blurred = median_filter(img_bin, size=100)
        img_final = img_bin - blurred.astype(float)
        fits.writeto(img_files[ii].replace('.fits', '_bin_nobkg.fits'), img_final, clobber=True)

    return
       
def find_stars_on_binned():
    # img_files = ['off_0S_0E_1.fits', 'off_0S_0E_2.fits', 'off_10S_0E_1.fits', 'off_10S_0E_2.fits',
    #              'off_10S_30E_1.fits', 'off_10S_30E_2.fits']
    img_files = ['off_0S_0E_1_bin_nobkg.fits']

    
    for ii in range(len(img_files)):
        print("Working on image: ", img_files[ii])
        img = fits.getdata(img_files[ii])

        # Calculate the bacgkround and noise (iteratively)
        print("\t Calculating background")
        threshold = 3
        for nn in range(5):
            if nn == 0:
                bkg_mean = img.mean()
                bkg_std = img.std()
            else:
                bkg_mean = img[good_pix].mean()
                bkg_std = img[good_pix].std()

            bad_hi = bkg_mean + (threshold * bkg_std)
            bad_lo = bkg_mean - (threshold * bkg_std)

            good_pix = np.where((img < bad_hi) & (img > bad_lo))
            print(nn, bkg_mean, bkg_std, threshold)
        
        bkg_mean = img[good_pix].mean()
        bkg_std = img[good_pix].std()
        print('Final', bkg_mean, bkg_std, threshold)
            
        threshold = (3.0 * bkg_std)

        print('\t Detecting Stars')
        daofind = DAOStarFinder(fwhm=4., threshold = threshold, exclude_border=True)
        sources = daofind(img - bkg_mean)

        # Save sources to output table
        pdb.set_trace()
        
    return sources




        

def find_outlier_pixels(data, tolerance=3, worry_about_edges=True, median_filter_size=10):
    """Find the hot or dead pixels in a 2D dataset. 

    Parameters
    ----------
    data: numpy array (2D)
        image
    tolerance: float (default 3)
        number of standard deviations used to cutoff the hot pixels
    worry_about_edges: boolean (default True)
        If you want to ignore the edges and greatly speed up the code, then set
        worry_about_edges to False.
        
    Return
    ------
    The function returns a list of hot pixels and also an image with with hot pixels removed

    """
    from scipy.ndimage import median_filter
    blurred = median_filter(data, size=median_filter_size)
    difference = data - blurred.astype(float)
    threshold = tolerance * np.std(difference)

    # Find the hot pixels, but ignore the edges
    hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1]) > threshold) )
    hot_pixels = np.array(hot_pixels) + 1 # because we ignored the first row and first column

    fixed_image = np.copy(data) # This is the image with the hot pixels removed
    for y, x in zip(hot_pixels[0], hot_pixels[1]):
        fixed_image[y, x] = blurred[y, x]

    if worry_about_edges == True:
        height, width = np.shape(data)

        ###Now get the pixels on the edges (but not the corners)###

        #left and right sides
        for index in range(1,height-1):
            #left side:
            med  = np.median(data[index-1:index+2,0:2])
            diff = np.abs(data[index,0] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[0]]  ))
                fixed_image[index,0] = med

            #right side:
            med  = np.median(data[index-1:index+2,-2:])
            diff = np.abs(data[index,-1] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[width-1]]  ))
                fixed_image[index,-1] = med

        #Then the top and bottom
        for index in range(1,width-1):
            #bottom:
            med  = np.median(data[0:2,index-1:index+2])
            diff = np.abs(data[0,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[0],[index]]  ))
                fixed_image[0,index] = med

            #top:
            med  = np.median(data[-2:,index-1:index+2])
            diff = np.abs(data[-1,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[height-1],[index]]  ))
                fixed_image[-1,index] = med

        ###Then the corners###

        #bottom left
        med  = np.median(data[0:2,0:2])
        diff = np.abs(data[0,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[0]]  ))
            fixed_image[0,0] = med

        #bottom right
        med  = np.median(data[0:2,-2:])
        diff = np.abs(data[0,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[width-1]]  ))
            fixed_image[0,-1] = med

        #top left
        med  = np.median(data[-2:,0:2])
        diff = np.abs(data[-1,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[0]]  ))
            fixed_image[-1,0] = med

        #top right
        med  = np.median(data[-2:,-2:])
        diff = np.abs(data[-1,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[width-1]]  ))
            fixed_image[-1,-1] = med

    return hot_pixels, fixed_image

    
        

        
