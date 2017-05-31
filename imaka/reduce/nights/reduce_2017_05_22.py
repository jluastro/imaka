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
import pdb

root_dir = '/Users/fatimaabdurrahman/Desktop/Research/RUN5/20170520/FLI/'

def make_flat():
    flat_raw_dir = root_dir + 'twilights/'
    dark_raw_dir = flat_raw_dir
    flat_out_dir = root_dir + "reduce/calib/"

    util.mkdir(flat_out_dir)
    
    flat_num = np.arange(106,  132)
    flat_frames = ['{0:s}twi_{1:04d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:04d}.fits'.format(dark_raw_dir, ss) for ss in flat_num]

    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat.fits')

    return

def make_sky():
    data_dir = root_dir + 'FLD2_2/'
    sky_dir = root_dir + 'reduce/sky/'
    
    util.mkdir(sky_dir)

    sky_num = np.arange(94, 106)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'FLD2_2_sky.fits')
    
    return
 
def reduce_FLD2():
    sky_dir = root_dir + 'reduce/sky/' 
    data_dir = root_dir + 'FLD2_2/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/FLD2_2/'
    
    util.mkdir(out_dir)

    
    # Open Loop
    fnum = [14, 17, 29, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59]
    fnum += [62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame=flat_dir + 'flat.fits')


    # Closed A
    fnum = [12, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55]
    fnum += [58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91]
    img_files = [data_dir + 'obj{0:04d}_cA.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed D
    fnum = [18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63]
    fnum += [66, 69, 72, 75, 78, 81, 84, 87, 90, 93]
    img_files = [data_dir + 'obj{0:04d}_cB.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    
    return
    
def find_stars_FLD2():
    data_dir = root_dir + 'reduce/FLD2_2/'

    # Open Loop
    fnum = [14, 17, 29, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59]
    fnum += [62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    #Closed A
    
    fnum = [12, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55]
    fnum += [58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91]
    img_files = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    # Closed D
    fnum = [18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63]
    fnum += [66, 69, 72, 75, 78, 81, 84, 87, 90, 93]
    img_files = [data_dir + 'obj{0:04d}_cD_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
    fnum = [14, 17, 29, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59]
    fnum += [62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed A
    fnum = [12, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55]
    fnum += [58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91]
    img_files = [data_dir + "obj{0:04d}_cA_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed D
    fnum = [18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63]
    fnum += [66, 69, 72, 75, 78, 81, 84, 87, 90, 93]
    img_files = [data_dir + "obj{0:04d}_cD_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedD.fits')
    
    return 

def stack_FLD2():

    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir + 'reduce/stats/'
    stacks_dir = root_dir + 'reduce/stacks/'

    util.mkdir(stacks_dir)

    # Open Loop
    open_img_nums = [14, 17, 29, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59]
    open_img_nums += [62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in open_img_nums]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in open_img_nums]
    open_output_root = stacks_dir + 'FLD2_2_stack_open'
    #reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed A
    closed_img_nums = [12, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55]
    closed_img_nums += [58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91]
    closed_images = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cA_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = stacks_dir + 'FLD2_2_stack_closedA'
    #reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed D
    closed_img_nums = [18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63]
    closed_img_nums += [66, 69, 72, 75, 78, 81, 84, 87, 90, 93]
    closed_images = [data_dir + 'obj{0:04d}_cD_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cD_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = stacks_dir + 'FLD2_2_stack_closedD'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return


def analyze_stacks():
    data_dir = root_dir + 'reduce/stacks/'
    stats_dir = root_dir + 'reduce/stats/'
    
    img_files = [data_dir + 'FLD2_2_stack_open.fits',
                 data_dir + 'FLD2_2_stack_closedA.fits',
                 data_dir + 'FLD2_2_stack_closedD.fits']
    
    #Find stars in image
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
