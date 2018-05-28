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


root_dir = '//Volumes/DATA5/imaka/20180527/FLI/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilight/'
massdimm_dir = root_dir + 'reduce/massdimm/'

    
fnum_o  = [42, 46, 50, 54, 58, 62, 66, 70]
fnum_threeWFS_LS = [41, 45, 49, 53, 57, 61, 65, 69] 
fnum_threeWFSLS_B2_c = [43, 47, 51, 55, 59, 63, 67, 71] 
fnum_threeWFSMean_B2_c  = [40, 44, 48, 52, 56, 60, 64, 68, 72] 

   
#def make_flat(): 
# No darks, so made flat by hand using darks from a different night and twilights
# from this night.  Not ideal, but best option.  'flat.fits' in reduce/calib


def make_sky():
    util.mkdir(sky_dir)

    sky_num = np.arange(73, 82+1)
    sky_frames = ['{0:s}sky_{1:04d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_30_sky.fits')
    
    return
 
def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - threeWFS_LS
    img_files = [data_dir + 'obj{0:04d}_threeWFS_LS.fits'.format(ii) for ii in fnum_threeWFS_LS]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSLS_B2_c
    img_files = [data_dir + 'obj{0:04d}_threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - threeWFSMean_B2_c 
    img_files = [data_dir + 'obj{0:04d}_threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c ]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_30_sky.fits', flat_frame =flat_dir+"flat.fits")


    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.find_stars(img_files, fwhm=9, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
    
    #Closed Loop - 3S WFS
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - 3L WFS
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=20, bottom_slice=0)
                              
    #Closed Loop - 4 WFS
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
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
    
    #Closed Loop - 4 WFS
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed_threeWFSMean_B2_c.fits')

    return


def calc_mof_stats():

    # Open Loop
    stats_file = stats_dir + 'stats_open.fits'
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - 3S
    stats_file = stats_dir + 'stats_closed_threeWFS_LS.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - 3L
    stats_file = stats_dir + 'stats_closed_threeWFSLS_B2_c.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop - 4
    stats_file = stats_dir + 'stats_closed_threeWFSMean_B2_c.fits'
    img_files = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    moffat.fit_moffat(img_files, stats_file)

    return


def append_massdimm():

    massdimm.fetch_data('20180527', massdimm_dir)
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

    # Open Loop - 30 s
    open_images = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - 3S WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_threeWFS_LS_clean.fits'.format(ii) for ii in fnum_threeWFS_LS]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFS_LS_clean_stars.txt'.format(ii) for ii in fnum_threeWFS_LS]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFS_LS_clean'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSLS_B2_c_clean.txt'.format(ii) for ii in fnum_threeWFSLS_B2_c]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSLS_B2_c'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - 3L WFS - 30 S
    closed_images = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_clean.fits'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    closed_starlists = [out_dir + 'obj{0:04d}_threeWFSMean_B2_c_cleann_stars.txt'.format(ii) for ii in fnum_threeWFSMean_B2_c]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFSMean_B2_c_clean'
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
    
