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
import poppy
import pdb
from astroscrappy import detect_cosmics

scale = 40.0 # mas/pixel

def subtract_dark(image_file, dark_file):
    # Dark subtract an image.
    img, i_hdr = fits.getdata(image_file, header=True)
    drk, d_hdr = fits.getdata(dark_file, header=True)

    # Rescale the dark to match the integration time... this isn't perfect.
    drk_time = d_hdr['EXPTIME']
    img_time = i_hdr['EXPTIME']

    drk_scale = img_time / drk_time

    img_ds = img - (drk * drk_scale)

    return img_ds

def fix_bad_pixels(img):
    crmask, cleanarr = detect_cosmics(img, sigclip=5, sigfrac=0.3, objlim=5.0, gain=1.0,
                                      readnoise=6.5, satlevel=65535.0, pssl=0, niter=4,
                                      sepmed=True, cleantype='meanmask', fsmode='median',
                                      psfmodel='gauss', psffwhm=2.5, psfsize=7,
                                      psfk=None, psfbeta=4.765, verbose=False)
                       
    return cleanarr


