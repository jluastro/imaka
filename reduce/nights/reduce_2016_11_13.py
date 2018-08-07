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
from astropy.modeling import models, fitting
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

    reduce_fli.clean_images(img_files)

    return

       
def find_stars_on_binned():
    # img_files = ['off_0S_0E_1.fits', 'off_0S_0E_2.fits', 'off_10S_0E_1.fits', 'off_10S_0E_2.fits',
    #              'off_10S_30E_1.fits', 'off_10S_30E_2.fits']
    img_files = ['off_0S_0E_1_bin_nobkg.fits']

    
    reduce_fli.find_stars_bin10(img_files)

    return



        

