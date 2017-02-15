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

root_dir = '/Users/fatimaabdurrahman/Desktop/20170214/FLI/'

def make_flat():
    flat_raw_dir = root_dir + 'twilights/'
    dark_raw_dir = root_dir + 'dark/'
    flat_out_dir = root_dir + "reduce/calib/"

    #util.mkdir(flat_out_dir)
    
    flat_num = np.arange(0, 28)
    flat_frames = ['{0:s}twi_{1:04d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:04d}.fits'.format(dark_raw_dir, ss) for ss in flat_num]
    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat.fits')

    return


def make_sky():
    data_dir = root_dir + 'Pleiades/'
    sky_dir = root_dir + 'reduce/sky/'

    sky_num = np.arange(78, 82)
    sky_frames = ['{0:s}sky__tt{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky.fits')
    
    return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    #util.mkdir(out_dir)

    # Open Loop
    fnum1 = [24, 32, 35, 38, 41, 44, 48, 52, 55, 58, 61, 64, 67, 70, 73, 76]
    fnum2 = [85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 116, 119]
    fnum = fnum1 + fnum2
    img_files = [data_dir + 'obj_o{0:04d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop
    fnum1 = [31, 34, 37, 40, 42, 46, 47, 50, 51, 54, 57, 60, 63, 66, 69, 72]
    fnum2 = [75, 86, 89, 92, 95, 98, 101, 104, 107, 110, 113, 118]
    fnum = fnum1 + fnum2
    img_files = [data_dir + 'obj_c{0:04d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # tt
    fnum1 = [33, 36, 39, 43, 45, 49, 53, 56, 59, 62, 65, 68, 71, 74, 83, 84]
    fnum2 = [87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 115, 117, 120]
    fnum = fnum1 + fnum2
    img_files = [data_dir + 'obj_tt{0:04d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')
    

    return
    
    
def find_stars_pleiades():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    # Open loop
    fnum1 = [24, 32, 35, 38, 41, 44, 48, 52, 55, 58, 61, 64, 67, 70, 73, 76]
    fnum2 = [85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 116, 119]
    fnum = fnum1 + fnum2
    img_files = [reduce_dir + 'obj_o{0:04d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    # Closed loop
    fnum1 = [31, 34, 37, 40, 42, 46, 47, 50, 51, 54, 57, 60, 63, 66, 69, 72]
    fnum2 = [75, 86, 89, 92, 95, 98, 101, 104, 107, 110, 113, 118]
    fnum = fnum1 + fnum2    
    img_files = [reduce_dir + 'obj_c{0:04d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=4, threshold=6)
    
    # tt
    fnum1 = [33, 36, 39, 43, 45, 49, 53, 56, 59, 62, 65, 68, 71, 74, 83, 84]
    fnum2 = [87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 115, 117, 120]
    fnum = fnum1 + fnum2    
    img_files = [reduce_dir + 'obj_tt{0:04d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)
    

    return


def calc_star_stats():
    reduce_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

#      # Open loop
#     fnum1 = [24, 32, 35, 38, 41, 44, 48, 52, 55, 58, 61, 64, 67, 70, 73, 76]
#     fnum2 = [85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 116, 119]
#     fnum = fnum1 + fnum2
#     img_files = [reduce_dir + 'obj_o{0:04d}_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')
    
    # Closed loop
    fnum1 = [31, 34, 37, 40, 42, 46, 47, 50, 51, 54, 57, 60, 63, 66, 69, 72]
    fnum2 = [75, 86, 89, 92, 95, 98, 101, 104, 107, 110, 113, 118]
    fnum = fnum1 + fnum2    
    img_files = [reduce_dir + 'obj_c{0:04d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')
    
    # tt
    fnum1 = [33, 36, 39, 43, 45, 49, 53, 56, 59, 62, 65, 68, 71, 74, 83, 84]
    fnum2 = [87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 115, 117, 120]
    fnum = fnum1 + fnum2    
    img_files = [reduce_dir + 'obj_tt{0:04d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_tt.fits')
    
    return
