import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
import reduce_fli
import calib
import util
import pdb
import os

root_dir = '/Users/fatimaabdurrahman/Desktop/20170215/FLI/'


def make_sky():
    data_dir = root_dir + 'Pleiades/'
    sky_dir = root_dir + 'reduce/sky/'

    #util.mkdir(out_dir)

    sky_num = np.arange(89, 110)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky.fits')
    
    return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    #util.mkdir(out_dir)

    # Open Loop
    fnum = [44, 47, 50, 51, 53, 54, 56, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop
    fnum1 = [40, 41, 43, 45, 46, 48, 49, 52, 55, 57, 58, 59, 61, 62]
    fnum2 = [64, 67, 70, 73, 76, 79, 82, 85, 88]
    fnum = fnum1 + fnum2
    img_files_1 = [data_dir + 'obj{0:04d}_c.fits'.format(ii) for ii in fnum]
    
    fnum = np.arange(35, 40)
    img_files_2 = [data_dir + "obj{0:04d}.fits".format(ii) for ii in fnum]
    
    img_files = img_files_1 + img_files_2
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # tt
    fnum = [65, 68, 71, 74, 77, 80, 83, 86]
    img_files = [data_dir + 'obj{0:04d}_tt.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')
    

    return
    
    
def find_stars_pleiades():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    # Open loop
#     fnum = [44, 47, 50, 51, 53, 54, 56, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87]    
#     img_files = [reduce_dir + 'obj{0:04}_o_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    # Closed loop
    fnum1 = [40, 41, 43, 45, 46, 48, 49, 52, 55, 57, 58, 59, 61, 62]
    fnum2 = [64, 67, 70, 73, 76, 79, 82, 85, 88]
    fnum = fnum1 + fnum2
    img_files_1 = [reduce_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    
    fnum = np.arange(35, 40)
    img_files_2 = [reduce_dir + "obj{0:04d}_clean.fits".format(ii) for ii in fnum]
    
    img_files = img_files_1 + img_files_2    
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)

    # tt
    fnum = [65, 68, 71, 74, 77, 80, 83, 86]  
    img_files = [reduce_dir + 'obj{0:04d}_tt_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)
    

    return


def calc_star_stats():
    reduce_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

     # Open loop
    fnum = [44, 47, 50, 51, 53, 54, 56, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87]
    img_files = [reduce_dir + 'obj{0:04}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')
    
    # Closed loop
    fnum1 = [40, 41, 43, 45, 46, 48, 49, 52, 55, 57, 58, 59, 61, 62]
    fnum2 = [64, 67, 70, 73, 76, 79, 82, 85, 88]
    fnum = fnum1 + fnum2
    img_files_1 = [reduce_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    
    fnum = np.arange(35, 40)
    img_files_2 = [reduce_dir + "obj{0:04d}_clean.fits".format(ii) for ii in fnum]
    
    img_files = img_files_1 + img_files_2
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')
    
    # tt
    fnum = [65, 68, 71, 74, 77, 80, 83, 86]  
    img_files = [reduce_dir + 'obj{0:04d}_tt_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_tt.fits')
    
    return
