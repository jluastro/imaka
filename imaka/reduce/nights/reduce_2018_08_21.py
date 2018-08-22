import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
import matplotlib
matplotlib.use('Agg')
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
from imaka.analysis import moffat
import os, shutil
import pdb
from imaka.reduce import massdimm
from imaka.reduce import reduce_STA
import matplotlib

################################################################

root_dir = '//Volumes/DATA5/imaka/20180821/sta/'
sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

################################################################

# Pre-shimming
fnum_LS_nrej7_1  = [31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71,\
                  75, 79, 83, 87, 91, 95, 99, 103, 107, 111]
fnum_o_1         = [32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, \
                  76, 80, 84, 88, 92, 96, 100, 104, 108, 112]
fnum_LS_bin2_1   = [33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, \
                  77, 81, 85, 89, 93, 97, 101, 105, 109, 113]
fnum_Mean_bin2_1 = [34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, \
                  78, 82, 86, 90, 94, 98, 102, 106, 110, 114]

# Post-shimming
fnum_LS_nrej7_2  = [115, 119, 123, 136, 140, 144, 148, 152, 156, 160]
fnum_o_2         = [116, 120, 124, 137, 141, 145, 149, 153, 157, 161]
fnum_LS_bin2_2   = [117, 121, 125, 138, 142, 146, 150, 154, 158, 162]
fnum_Mean_bin2_2 = [118, 122, 126, 139, 143, 147, 151, 155, 159, 163]

################################################################

def make_flat(): 

    util.mkdir(flat_dir)
    flat_num = np.arange(0, 9+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    calib.makeflat(scan_flat_frames, [], flat_dir + 'flat.fits', darks=False)

    return


def make_sky():

    util.mkdir(sky_dir)
    sky_num = np.arange(164, 171+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'FLD2_sky.fits')
    
    return


def reduce_FLD2():

    util.mkdir(out_dir)

    #Pre-shimming:
    
    # Closed - LS_nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c.fits'.format(ii) for ii in fnum_LS_nrej7_1]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_LS_nrej7_1]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Open 
    img_files = [data_dir + 'obj{0:03d}_o.fits'.format(ii) for ii in fnum_o_1]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o_1]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir+"flat.fits")
                
    # Closed - LS_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_LS_bin2_1]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan.fits'.format(ii) for ii in fnum_LS_bin2_1]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")
                
    # Closed - Mean_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_Mean_bin2_1]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan.fits'.format(ii) for ii in fnum_Mean_bin2_1]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    #Post-shimming:
    
    # Closed - LS_nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c.fits'.format(ii) for ii in fnum_LS_nrej7_2]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_LS_nrej7_2]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Open 
    img_files = [data_dir + 'obj{0:03d}_o.fits'.format(ii) for ii in fnum_o_2]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o_2]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir+"flat.fits")
                
    # Closed - LS_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_LS_bin2_2]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan.fits'.format(ii) for ii in fnum_LS_bin2_2]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")
                
    # Closed - Mean_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_Mean_bin2_2]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan.fits'.format(ii) for ii in fnum_Mean_bin2_2]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")
                
    return



def find_stars_FLD2():

    # Pre-shimming

    # Closed - LS_nrej7
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_LS_nrej7_1]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Open 
    img_files = [data_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o_1]
    reduce_fli.find_stars(img_files, fwhm=10, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Closed - LS_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_LS_bin2_1]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Closed - Mean_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan_clean.fits'.format(ii) for ii in fnum_Mean_bin2_1]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Post-shimming

    # Closed - LS_nrej7
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_LS_nrej7_2]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Open 
    img_files = [data_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o_2]
    reduce_fli.find_stars(img_files, fwhm=10, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Closed - LS_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_LS_bin2_2]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    # Closed - Mean_bin2
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan_clean.fits'.format(ii) for ii in fnum_Mean_bin2_2]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                               left_slice =20, right_slice=20, top_slice=30, bottom_slice=30)

    return


def calc_star_stats():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    stats_file = stats_dir + 'stats_open.fits'
    #reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    #moffat.fit_moffat(img_files, stats_file)

    # EXPERIMENT 1:
    
    # Closed - LS_c
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_LS_c]
    stats_file = stats_dir + 'stats_closed_LS.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    # Closed - LS_B2_c
    img_files = [out_dir + 'obj{0:03d}threeWFSLS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_LS_B2_c]
    stats_file = stats_dir + 'stats_closed_LS_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - Mean_B2_c
    img_files = [out_dir + 'obj{0:03d}threeWFSMean_B2_c_scan_clean.fits'.format(ii) for ii in fnum_Mean_B2_c]
    stats_file = stats_dir + 'stats_closed_Mean_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # EXPERIMENT 2

    # Closed - nrej1
    img_files = [out_dir + 'obj{0:03d}threeWFSLS_nrej1_c_scan_clean.fits'.format(ii) for ii in fnum_nrej1]
    stats_file = stats_dir + 'stats_closed_nrej1.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - nrej4
    img_files = [out_dir + 'obj{0:03d}threeWFSLS_nrej4_c_scan_clean.fits'.format(ii) for ii in fnum_nrej4]
    stats_file = stats_dir + 'stats_closed_nrej4.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - nrej7
    img_files = [out_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_nrej7]
    stats_file = stats_dir + 'stats_closed_nrej7.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    return


def append_massdimm():

    massdimm.fetch_data('20180602', massdimm_dir)
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

    # EXPERIMENT 1

    # Closed - LS_c
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_LS_c]
    starlists = [out_dir  + 'obj{0:03d}threeWFS_LS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_LS_c]
    output_root = stacks_dir + 'FLD2_stack_closed_LS'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

    # Closed - LS_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_LS_B2_c]
    starlists = [out_dir  + 'obj{0:03d}threeWFSLS_B2_c_scan_clean_stars.txt'.format(ii) for ii in fnum_LS_c]
    output_root = stacks_dir + 'FLD2_stack_closed_LS_B2'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

    # Closed - Mean_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan_clean.fits'.format(ii) for ii in fnum_Mean_B2_c]
    starlists = [out_dir  + 'obj{0:03d}threeWFSMean_B2_c_scan_clean_stars.txt'.format(ii) for ii in fnum_Mean_B2_c]
    output_root = stacks_dir + 'FLD2_stack_closed_Mean_B2'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

    # EXPERIMENT 2

    # Closed - nrej1
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej1_c_scan_clean.fits'.format(ii) for ii in fnum_nrej1]
    starlists = [out_dir  + 'obj{0:03d}threeWFSLS_nrej1_c_scan_clean_stars.txt'.format(ii) for ii in fnum_nrej1]
    output_root = stacks_dir + 'FLD2_stack_closed_nrej1'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

    # Closed - nrej4
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej4_c_scan_clean.fits'.format(ii) for ii in fnum_nrej4]
    starlists = [out_dir  + 'obj{0:03d}threeWFSLS_nrej4_c_scan_clean_stars.txt'.format(ii) for ii in fnum_nrej4]
    output_root = stacks_dir + 'FLD2_stack_closed_nrej4'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

    # Closed - nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan_clean.fits'.format(ii) for ii in fnum_nrej7]
    starlists = [out_dir  + 'obj{0:03d}threeWFS_LS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_nrej7]
    output_root = stacks_dir + 'FLD2_stack_closed_nrej7'
    reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')

def analyze_stacks():

    open_img_files = [stacks_dir + 'FLD2_stack_open.fits']

    closed_img_files = [stacks_dir + 'FLD2_stack_closed_LS.fits', \
                        stacks_dir + 'FLD2_stack_closed_LS_B2.fits', \
                        stacks_dir + 'FLD2_stack_closed_Mean_B2.fits', \
                        stacks_dir + 'FLD2_stack_closed_nrej1.fits', \
                        stacks_dir + 'FLD2_stack_closed_nrej4.fits', \
                        stacks_dir + 'FLD2_stack_closed_nrej7.fits']
    
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
    
