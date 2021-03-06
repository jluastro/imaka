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
from imaka.reduce import reduce_STA
import matplotlib
matplotlib.use('Agg')

root_dir = '//Volumes/DATA5/imaka/20180601/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

fnum_o = [23, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, \
              58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, \
              90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116,\
              118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 144,\
              146, 148, 150, 152, 154, 156, 158, 160]
fnum_c = [22, 24, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, \
              55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, \
              87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 111, 113, 115,  \
              117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141,\
              143, 145, 147, 149, 151, 153, 155, 157, 159]


def make_flat(): 

    util.mkdir(flat_dir)
    
    flat_num = np.arange(171,  182+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    calib.makedark(scan_flat_frames, flat_dir + 'flat.fits')

    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(161, 170+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'FLD2_sky.fits')
    
    return


def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'obj{0:03d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed - threeWFS_LS
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c.fits'.format(ii) for ii in fnum_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")


    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.find_stars(img_files, fwhm=10, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    
    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    stats_file = stats_dir + 'stats_open.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    stats_file = stats_dir + 'stats_closed_threeWFS_LS_c.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    return


def append_massdimm():

    massdimm.fetch_data('20180601', massdimm_dir)
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
    open_images = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - threeWFS_LS
    closed_images = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    closed_starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c]
    closed_output_root = stacks_dir + 'FLD2_stack_threeWFS_LS_c'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    return


def analyze_stacks():

    open_img_files = [stacks_dir + 'FLD2_stack_open.fits']

    closed_img_files = [stacks_dir + 'FLD2_stack_threeWFS_LS_c.fits']
    
    #Find stars in image
    reduce_fli.find_stars(open_img_files, fwhm=10, threshold=100, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    reduce_fli.find_stars(closed_img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
        
    # Calc stats on all the stacked images
    #reduce_fli.calc_star_stats(open_img_files+closed_img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
def split_filt():

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.fourfilt(img_files, starlists)

    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.fourfilt(img_files, starlists)

    return

def calc_star_stats_fourfilt():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    

    # R
    stats_file = stats_dir + 'stats_open_R.fits'
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_R_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='R')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # V
    stats_file = stats_dir + 'stats_open_V.fits'
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_V_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='V')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # B
    stats_file = stats_dir + 'stats_open_B.fits'
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_B_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='B')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # I
    stats_file = stats_dir + 'stats_open_I.fits'
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_I_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='I')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_c]

    # R
    stats_file = stats_dir + 'stats_closed_R.fits'
    starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_R_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='R')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # V
    stats_file = stats_dir + 'stats_closed_V.fits'
    starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_V_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='V')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # B
    stats_file = stats_dir + 'stats_closed_B.fits'
    starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_B_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='B')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # I
    stats_file = stats_dir + 'stats_closed_I.fits'
    starlists = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean_I_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='I')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    return
