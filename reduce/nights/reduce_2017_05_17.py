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
import os, shutil


root_dir = '/Users/fatimaabdurrahman/Desktop/Research/RUN5/20170517/FLI/'


def make_flat():
    flat_raw_dir = root_dir + 'twilights/'
    dark_raw_dir = root_dir + 'dark/'
    flat_out_dir = root_dir + "reduce/calib/"

    util.mkdir(flat_out_dir)
    
    flat_num = np.arange(125, 152)
    flat_frames = ['{0:s}twi_{1:04d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:04d}.fits'.format(dark_raw_dir, ss) for ss in flat_num]
    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat.fits')

    return


def make_sky():
    data_dir = root_dir + 'FLD2/'
    sky_dir = root_dir + 'reduce/sky/'
    
    util.mkdir(sky_dir)

    sky_num = np.arange(116, 124)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_sky.fits')
    
    return


def reduce_FLD2():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'FLD2/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/FLD2/'
    
    util.mkdir(out_dir)

    # Open Loop
    fnum = [63, 67, 71, 75, 80, 84, 88, 92, 96,  100, 104, 108, 112, 124]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop
    fnum = [64, 68, 72, 76, 77, 81, 85, 89, 93,97, 101, 105, 109, 113]
    img_files = [data_dir + 'obj{0:04d}_c.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed A
    fnum = [65, 69, 73, 78, 82, 86,90, 94,98, 102, 106, 110,114]
    img_files = [data_dir + 'obj{0:04d}_cA.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    # Closed B
    fnum = [66, 70, 74,79, 83, 87, 91, 95, 99, 103, 107, 111,115]
    img_files = [data_dir + 'obj{0:04d}_cB.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    return


def find_stars_FLD2():
    data_dir = root_dir + 'reduce/FLD2/'

    # Open Loop
    fnum = [63, 67, 71, 75, 80, 84, 88, 92, 96,  100, 104, 108, 112, 124]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    # Closed Loop
    
    fnum = [64, 68, 72, 76, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)
    
    #Closed A
    
    fnum = [65, 69, 73, 78, 82, 86,90, 94,98, 102, 106, 110,114]
    img_files = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    
    
    # Closed B
    
    fnum = [66, 70, 74,79, 83, 87, 91, 95, 99, 103, 107, 111,115]
    img_files = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/FLD2/'
    stats_dir = root_dir +'reduce/stats/'

    
     # Open Loop
    fnum = [63, 67, 71, 75, 80, 84, 88, 92, 96, 100, 104, 108, 112, 124]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop
    fnum = [64, 68, 72, 76, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')

    # Closed A
    fnum = [65, 69, 73, 78, 82, 86,90, 94, 98, 102, 106, 110,114]
    img_files = [data_dir + "obj{0:04d}_cA.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed B
    fnum = [66, 70, 74,79, 83, 87, 91, 95, 99, 103, 107, 111,115]
    img_files = [data_dir + "obj{0:04d}_cB.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedB.fits')
    
    return

def stack_FLD2():

    data_dir = root_dir + 'reduce/FLD2/'
    stats_dir = root_dir +'reduce/stats/'

    # Open Loop
    open_img_nums = [63, 67, 71, 75, 80, 84, 88, 92, 96, 100, 104, 108, 112, 124]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in open_img_nums]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in open_img_nums]
    open_output_root = data_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop
    closed_img_nums = [64, 68, 72, 76, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113]
    closed_images = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_c_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = data_dir + 'FLD2_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed A
    closed_img_nums = [65, 69, 73, 78, 82, 86, 90, 94,98, 102, 106, 110,114]
    closed_images = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cA_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = data_dir + 'FLD2_stack_closedA'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed B
    closed_img_nums = [66, 70, 74,79, 83, 87, 91, 95, 99, 103, 107, 111,115]
    closed_images = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cB_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = data_dir + 'FLD2_stack_closedB'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return


def analyze_stacks():
    data_dir = root_dir + 'reduce/stacks/'
    
    img_files = [data_dir + 'FLD2_stack_open.fits',
                 data_dir + 'FLD2_stack_closed.fits',
                 data_dir + 'FLD2_stack_closedA.fits',
                 data_dir + 'FLD2_stack_closedB.fits']
    
    #Find stars in image
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats='stats_stacks.fits')

    return
    
