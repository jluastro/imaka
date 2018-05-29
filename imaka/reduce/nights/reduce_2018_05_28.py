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
from imaka.reduce import massdimm


root_dir = '//Volumes/DATA5/imaka/20180528/FLI/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilight/'
massdimm_dir = root_dir + 'reduce/massdimm/'

fnum_o  = [111, 115, 119, 123, 127, 131, 135, 139, 143, 147, 151]
fnum_threeWFS_LS = [109, 110, 114, 118, 122, 126, 130, 134, 138, 142, 146, 150] 
fnum_threeWFSLS_B2_c = [112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152] 
fnum_threeWFSMean_B2_c  = [113, 117, 121, 125, 129, 133, 137, 141, 145, 149, 153]

fnum_o_2filt = [166, 170, 174, 178, 182, 186]
fnum_threeWFS_LS_2filt = [165, 169, 173, 177, 181, 185]
fnum_threeWFSLS_B2_c_2filt = [167, 171, 175, 179, 183, 187]
fnum_threeWFSMean_B2_c_2filt = [164, 168, 172, 176, 180, 184, 188]

   
#def make_flat(): 
#Using flat from 20180527, copied from
#/Volumes/DATA5/imaka/20180527/FLI/reduce/calib


def make_sky():
    util.mkdir(sky_dir)

    sky_num = np.arange(154, 163+1)
    sky_frames = ['{0:s}sky_{1:04d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_sky.fits')

    sky_num = np.arange(189, 198+1)
    sky_frames = ['{0:s}sky_{1:04d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_2filt_sky.fits')
    
    return
 
def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - threeWFS_LS
    img_files = [data_dir + 'obj{0:04d}_threeWFS_LS.fits'.format(ii) for ii in fnum_threeWFS_LS]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSLS_B2_c
    img_files = [data_dir + 'obj{0:04d}_threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSMean_B2_c 
    img_files = [data_dir + 'obj{0:04d}_threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c ]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Open Loop - 2 filter
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o_2filt]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2filt_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - threeWFS_LS - 2 filter
    img_files = [data_dir + 'obj{0:04d}_threeWFS_LS.fits'.format(ii) for ii in fnum_threeWFS_LS_2filt]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2filt_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSLS_B2_c - 2 filter
    img_files = [data_dir + 'obj{0:04d}_threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c_2filt]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2filt_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSMean_B2_c - 2 filter
    img_files = [data_dir + 'obj{0:04d}_threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c_2filt]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2filt_sky.fits', flat_frame =flat_dir+"flat.fits")


    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o+fnum_o_2filt]
    reduce_fli.find_stars(img_files, fwhm=9, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
    
    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS+fnum_threeWFS_LS_2filt]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - threeWFSLS_B2_c
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c+fnum_threeWFSLS_B2_c_2filt]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - threeWFSMean_B2_c
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c+fnum_threeWFSMean_B2_c_2filt]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)

    return


def calc_star_stats():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFS_LS.fits')

    #Closed Loop - threeWFSLS_B2_c
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFSLS_B2_c.fits')
    
    #Closed Loop - threeWFSMean_B2_c
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFSMean_B2_c.fits')

    # Open Loop - 2 filter
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_2filt]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open_2filt.fits')

    #Closed Loop - threeWFS_LS - 2 filter
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS_2filt]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFS_LS_2filt.fits')

    #Closed Loop - threeWFSLS_B2_c - 2 filter
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c_2filt]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFSLS_B2_c_2filt.fits')
    
    #Closed Loop - threeWFSMean_B2_c - 2 filter
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c_2filt]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFSMean_B2_c_2filt.fits')

    return


def calc_mof_stats():

    # Open Loop
    stats_file = stats_dir + 'stats_open.fits'
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - threeWFS_LS
    stats_file = stats_dir + 'stats_closed_threeWFS_LS.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - threeWFSLS_B2_c
    stats_file = stats_dir + 'stats_closed_threeWFSLS_B2_c.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - threeWFSMean_B2_c
    stats_file = stats_dir + 'stats_closed_threeWFSMean_B2_c.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    moffat.fit_moffat(img_files, stats_file)

    # Open Loop - 2 filter
    stats_file = stats_dir + 'stats_open_2filt.fits'
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_2filt]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - threeWFS_LS  - 2 filter
    stats_file = stats_dir + 'stats_closed_threeWFS_LS_2filt.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS_2filt]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - 3LthreeWFSLS_B2_c - 2 filter
    stats_file = stats_dir + 'stats_closed_threeWFSLS_B2_c_2filt.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c_2filt]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - threeWFSMean_B2_c - 2 filter
    stats_file = stats_dir + 'stats_closed_threeWFSMean_B2_c_2filt.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c_2filt]
    moffat.fit_moffat(img_files, stats_file)

    return


def append_massdimm():

    massdimm.fetch_data('20180528', massdimm_dir)
    stats_tables = glob.glob(root_dir + 'reduce/stats/stats*.fits')

    for stats in stats_tables:
        if 'mdp.fits' not in stats:
            print('Adding MASS/DIMM to ' + stats)
            massdimm.append_mass_dimm(stats, massdimm_dir)
        else:
            print('Skipping ' + stats)

    return


def stack_FLD2():

    util.mkdir(stacks_dir)

    # Open Loop
    open_images = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - threeWFS_LS
    closed_images = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFS_LS_clean_stars.txt'.format(ii) for ii in fnum_threeWFS_LS]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFS_LS'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - threeWFSLS_B2_c
    closed_images = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean_stars.txt'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSLS_B2_c'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - threeWFSMean_B2_c
    closed_images = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean_stars.txt'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSMean_B2_c'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Open Loop - 2 filter
    open_images = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o_2filt]
    open_starlists = [out_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o_2filt]
    open_output_root = stacks_dir + 'FLD2_stack_open_2filt'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - threeWFS_LS - 2 filter
    closed_images = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS_2filt]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFS_LS_clean_stars.txt'.format(ii) for ii in fnum_threeWFS_LS_2filt]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFS_LS_2filt'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - threeWFSLS_B2_c - 2 filter
    closed_images = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c_2filt]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean_stars.txt'.format(ii) for ii in fnum_threeWFSLS_B2_c_2filt]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSLS_B2_c_2filt'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - threeWFSMean_B2_c - 2 filter
    closed_images = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c_2filt]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean_stars.txt'.format(ii) for ii in fnum_threeWFSMean_B2_c_2filt]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSMean_B2_c_2filt'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    
    return


def analyze_stacks():

    open_img_files = [stacks_dir + 'FLD2_stack_open.fits',
                      stacks_dir + 'FLD2_stack_open_2filt.fits']

    closed_img_files = [stacks_dir + 'FLD2_stack_threeWFSLS_B2_c.fits',
                        stacks_dir + 'FLD2_stack_threeWFSMean_B2_c.fits',
                        stacks_dir + 'FLD2_stack_threeWFS_LS.fits',
                        stacks_dir + 'FLD2_stack_threeWFSLS_B2_c_2filt.fits',
                        stacks_dir + 'FLD2_stack_threeWFSMean_B2_c_2filt.fits',
                        stacks_dir + 'FLD2_stack_threeWFS_LS_2filt.fits']
    
    #Find stars in image
    reduce_fli.find_stars(open_img_files, fwhm=9, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
    reduce_fli.find_stars(closed_img_files, fwhm=4, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
        
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(open_img_files+closed_img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
