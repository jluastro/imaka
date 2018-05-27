import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
from imaka.analysis import moffat
import os, shutil
import pdb

root_dir = '//Volumes/DATA5/imaka/20180526/FLI/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilight/'
    
fnum_o_30  = [5, 8, 11, 14, 18, 22, 26, 30, 34, 38, 42, 46]
fnum_c_3S_30 = [6, 9, 12, 15, 20, 24, 28, 32, 36, 40, 44, 48] # Closed loop, 3 WFS small
fnum_c_3L_30 = [16, 19, 23, 27, 31, 35, 39, 43, 47] # Closed loop, 3 WFS large
fnum_c_4_30  = [3, 4, 7, 10, 13, 17, 21, 25, 29, 33, 37, 41, 45] # Closed loop, 4 WFS
fnum_o_60    = [50, 54, 58, 65]
fnum_c_3S_60 = [52, 56, 60, 64]
fnum_c_3L_60 = [51, 55, 59, 63]
fnum_c_4_60  = [49, 53, 57, 61, 62] 
   
#def make_flat():
#    util.mkdir(flat_dir)
#    
#    flat_num = np.arange(10,  24+1)
#    flat_frames = ['{0:s}twi{1:02d}.fits'.format(twi_dir, ss) for ss in flat_num]
#    dark_frames = ['{0:s}dark{1:02d}.fits'.format(twi_dir, ss) for ss in flat_num]

#    calib.makeflat(flat_frames, dark_frames, flat_dir + 'flat.fits')

#    return


def make_sky():
    util.mkdir(sky_dir)

    sky_num = np.arange(66, 69+1)
    sky_frames = ['{0:s}sky_{1:04d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_30_sky.fits')

    sky_num = np.arange(70, 75+1)
    sky_frames = ['{0:s}sky_{1:04d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_60_sky.fits')
    
    return
 
def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop - 30s
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o_30]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - 3 WFS Small - 30s
    img_files = [data_dir + 'obj{0:04d}_threewfs_small_c.fits'.format(ii) for ii in fnum_c_3S_30]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - 3 WFS Large - 30s
    img_files = [data_dir + 'obj{0:04d}_threeWFS_big_c.fits'.format(ii) for ii in fnum_c_3L_30]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - 4 WFS - 30s
    img_files = [data_dir + 'obj{0:04d}_fourWFS_c.fits'.format(ii) for ii in fnum_c_4_30]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")

     # Open Loop - 60s
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o_60]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_60_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - 3 WFS Small - 60s
    img_files = [data_dir + 'obj{0:04d}_threewfs_small_c.fits'.format(ii) for ii in fnum_c_3S_60]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_60_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - 3 WFS Large - 60s
    img_files = [data_dir + 'obj{0:04d}_threeWFS_big_c.fits'.format(ii) for ii in fnum_c_3L_60]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_60_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - 4 WFS - 60s
    img_files = [data_dir + 'obj{0:04d}_fourWFS_c.fits'.format(ii) for ii in fnum_c_4_60]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_60_sky.fits', flat_frame =flat_dir+"flat.fits")

    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_30+fnum_o_60]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
    
    #Closed Loop - 3S WFS
    img_files = [out_dir + 'obj{0:04d}_threewfs_small_c_clean.fits'.format(ii) for ii in fnum_c_3S_30+fnum_c_3S_60]
    reduce_fli.find_stars(img_files, fwhm=4, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - 3L WFS
    img_files = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean.fits'.format(ii) for ii in fnum_c_3L_30_fnum_c_3L_60]
    reduce_fli.find_stars(img_files, fwhm=4, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - 4 WFS
    img_files = [out_dir + 'obj{0:04d}_fourWFS_c_clean.fits'.format(ii) for ii in fnum_c_4_30+fnum_c_4_60]
    reduce_fli.find_stars(img_files, fwhm=4, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    return


def calc_star_stats():
    
    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_30+fnum_o_60]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    #Closed Loop - 3S WFS
    img_files = [out_dir + 'obj{0:04d}_threewfs_small_c_clean.fits'.format(ii) for ii in fnum_c_3S_30+fnum_c_3S_60]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_3WFS_S.fits')

    #Closed Loop - 3L WFS
    #img_files = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean.fits'.format(ii) for ii in fnum_c_3L_30_fnum_c_3L_60]
    #reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_3WFS_L.fits')
    
    #Closed Loop - 4 WFS
    #img_files = [out_dir + 'obj{0:04d}_fourWFS_c_clean.fits'.format(ii) for ii in fnum_c_4_30+fnum_c_4_60]
    #reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_4WFS.fits')

    return


def calc_mof_stats():

    # Open Loop
    stats_file = stats_dir + 'stats_open.fits'
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop
    stats_file = stats_dir + 'stats_closed.fits'
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum_c]
    moffat.fit_moffat(img_files, stats_file)


def stack_FLD2():

    util.mkdir(stacks_dir)

    # Open Loop - 30 s
    open_images = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_30]
    open_starlists = [out_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o_30]
    open_output_root = stacks_dir + 'FLD2_stack_open_30'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - 3S WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_threewfs_small_c_clean.fits'.format(ii) for ii in fnum_c_3S_30]
    closed_starlists = [out_dir + 'obj{0:04d}_threewfs_small_c_clean_stars.txt'.format(ii) for ii in fnum_c_3S_30]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_3S_30'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean.fits'.format(ii) for ii in fnum_c_3L_30]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean_stars.txt'.format(ii) for ii in fnum_c_3L_30]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_3L_30'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_fourWFS_c_clean.fits'.format(ii) for ii in fnum_c_4_30]
    closed_starlists = [out_dir + 'obj{0:04d}_fourWFS_c_clean_stars.txt'.format(ii) for ii in fnum_c_4_30]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_4_30'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Open Loop - 60 s
    open_images = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_60]
    open_starlists = [out_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o_60]
    open_output_root = stacks_dir + 'FLD2_stack_open_60'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - 3S WFS - 60 S
    closed_images = [out_dir + 'obj{0:04d}_threewfs_small_c_clean.fits'.format(ii) for ii in fnum_c_3S_60]
    closed_starlists = [out_dir + 'obj{0:04d}_threewfs_small_c_clean_stars.txt'.format(ii) for ii in fnum_c_3S_60]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_3S_60'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 60 S
    closed_images = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean.fits'.format(ii) for ii in fnum_c_3L_60]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFS_big_c_clean_stars.txt'.format(ii) for ii in fnum_c_3L_60]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_3L_60'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 60 S
    closed_images = [out_dir + 'obj{0:04d}_fourWFS_c_clean.fits'.format(ii) for ii in fnum_c_4_60]
    closed_starlists = [out_dir + 'obj{0:04d}_fourWFS_c_clean_stars.txt'.format(ii) for ii in fnum_c_4_60]
    closed_output_root = stacks_dir + 'FLD2_stack_closed_4_60'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    return


def analyze_stacks():

    open_img_files = [stacks_dir + 'FLD2_stack_open_30.fits',
                     stacks_dir + 'FLD2_stack_open_60.fits']

    closed_img_files = [stacks_dir + 'FLD2_stack_closed_3S_30.fits',
                       stacks_dir + 'FLD2_stack_closed_3L_30.fits',
                       stacks_dir + 'FLD2_stack_closed_4_30.fits',
                       stacks_dir + 'FLD2_stack_closed_3S_60.fits',
                       stacks_dir + 'FLD2_stack_closed_3L_60.fits',
                       stacks_dir + 'FLD2_stack_closed_4_60.fits']
    
    #Find stars in image
    reduce_fli.find_stars(open_img_files, fwhm=8, threshold=10, N_passes=2, plot_psf_compare=False)
    reduce_fli.find_stars(closed_img_files, fwhm=4, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(open_img_files+closed_img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
