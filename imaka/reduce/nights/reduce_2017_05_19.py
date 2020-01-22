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

root_dir = '/Users/fatimaabdurrahman/Desktop/Research/RUN5/20170519/FLI/'

#flat taken from 20170518, no twilights taken this night

def make_sky():
    data_dir = root_dir + 'FLD2_2/'
    sky_dir = root_dir + 'reduce/sky/'
    
    util.mkdir(sky_dir)

    sky_num = np.arange(172, 185)
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
    fnum = [4, 5, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45]
    fnum += [51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 85, 88, 89]
    fnum += [92, 95, 99, 103, 107, 110, 113, 116, 119, 122, 125, 128]
    fnum += [131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop
    fnum = [96, 100, 104]
    img_files = [data_dir + 'obj{0:04d}_c.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed A
    fnum = [6, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53]
    fnum += [56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 87, 91, 94, 98, 102, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145]
    fnum += [148, 151, 154, 157, 160, 163, 166, 169]
    img_files = [data_dir + 'obj{0:04d}_cA.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed B
    fnum = [10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58]
    fnum += [61, 64, 67, 70, 73, 76, 79, 82, 86, 90, 93, 97, 101, 105, 108, 111]
    fnum += [114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153]
    fnum += [156, 159, 162, 165, 168, 171]
    img_files = [data_dir + 'obj{0:04d}_cB.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    
    return
    
def find_stars_FLD2():
    data_dir = root_dir + 'reduce/FLD2_2/'

    # Open Loop
    #fnum = [4, 5, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45]
    #fnum += [51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 85, 88, 89]
    #fnum += [92, 95, 99, 103, 107, 110, 113, 116, 119, 122, 125, 128]
    fnum = [125, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    # Closed Loop
    
    fnum = [96, 100, 104]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)
    
    #Closed A
    
    fnum = [6, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53]
    fnum += [56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 87, 91, 94, 98, 102, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145]
    fnum += [148, 151, 154, 157, 160, 163, 166, 169]
    img_files = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    # Closed B
    fnum = [10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58]
    fnum += [61, 64, 67, 70, 73, 76, 79, 82, 86, 90, 93, 97, 101, 105, 108, 111]
    fnum += [114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153]
    fnum += [156, 159, 162, 165, 168, 171]
    img_files = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir +'reduce/stats/'

    util.mkdir(stats_dir)
    
    # Open Loop
    fnum = [4, 5, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45]
    fnum += [51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 85, 88, 89]
    fnum += [92, 95, 99, 103, 107, 110, 113, 116, 119, 122, 125, 128]
    fnum += [131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop
    fnum = [96, 100, 104]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')

    # Closed A
    fnum = [6, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53]
    fnum += [56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 87, 91, 94, 98, 102, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145]
    fnum += [148, 151, 154, 157, 160, 163, 166, 169]
    img_files = [data_dir + "obj{0:04d}_cA_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed B
    fnum = [10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58]
    fnum += [61, 64, 67, 70, 73, 76, 79, 82, 86, 90, 93, 97, 101, 105, 108, 111]
    fnum += [114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153]
    fnum += [156, 159, 162, 165, 168, 171]
    img_files = [data_dir + "obj{0:04d}_cB_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedB.fits')
    
    return 

def stack_FLD2():

    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir + 'reduce/stats/'
    stacks_dir = root_dir + 'reduce/stacks/'

    util.mkdir(stacks_dir)

    # Open Loop
    open_img_nums = [4, 5, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45]
    open_img_nums += [51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 85, 88, 89]
    open_img_nums += [92, 95, 99, 103, 107, 110, 113, 116, 119, 122, 125, 128]
    open_img_nums += [131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in open_img_nums]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in open_img_nums]
    open_output_root = stacks_dir + 'FLD2_2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop
    closed_img_nums = [96, 100, 104]
    closed_images = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_c_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = stacks_dir + 'FLD2_2_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed A
    closed_img_nums = [6, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53]
    closed_img_nums += [56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 87, 91, 94, 98, 102, 106]
    closed_img_nums += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145]
    closed_img_nums += [148, 151, 154, 157, 160, 163, 166, 169]
    closed_images = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cA_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = stacks_dir + 'FLD2_2_stack_closedA'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed B
    closed_img_nums = [10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58]
    closed_img_nums += [61, 64, 67, 70, 73, 76, 79, 82, 86, 90, 93, 97, 101, 105, 108, 111]
    closed_img_nums += [114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153]
    closed_img_nums += [156, 159, 162, 165, 168, 171]
    closed_images = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_cB_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = stacks_dir + 'FLD2_2_stack_closedB'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return


def analyze_stacks():
    data_dir = root_dir + 'reduce/stacks/'
    stats_dir = root_dir + 'reduce/stats/'
    
    img_files = [data_dir + 'FLD2_2_stack_open.fits',
                 data_dir + 'FLD2_2_stack_closed.fits',
                 data_dir + 'FLD2_2_stack_closedA.fits',
                 data_dir + 'FLD2_2_stack_closedB.fits']
    
    #Find stars in image
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
