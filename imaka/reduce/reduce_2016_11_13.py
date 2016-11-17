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
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import Background, detect_threshold, detect_sources
from photutils import source_properties, properties_table
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


    img_files = ['off_0S_0E_1.fits', 'off_0S_0E_2.fits', 'off_10S_0E_1.fits', 'off_10S_0E_2.fits',
                 'off_10S_30E_1.fits', 'off_10S_30E_2.fits']


    for ii in range(len(img_files)):
        img_ds = reduce_fli.dark_subtract(img_files[ii], dark_file)

        

    
        
