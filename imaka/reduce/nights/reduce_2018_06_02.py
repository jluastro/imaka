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

################################################################

root_dir = '//Volumes/DATA5/imaka/20180602/sta/'
sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

################################################################

# Open Loop
fnum_o         = [87, 91, 95, 99, 103, 107, 111, 115, 119, \
                  123, 127, 132, 132, 136, 140, 150, 154, 158]

# Experiment 1
fnum_LS_c      = [83, 84, 85, 86, 90, 94, 98, 102, 106, 110, 114]
fnum_LS_B2_c   = [88, 92, 96, 100, 104, 108, 112, 116]
fnum_Mean_B2_c = [89, 93, 97, 101, 105, 109, 113, 117]

#Experiment 2

fnum_nrej1 = [120, 124, 128, 133, 137, 141, 151, 155, 159]
fnum_nrej4 = [121, 125, 134, 138, 142, 152, 156, 160]
fnum_nrej7 = [118, 122, 126, 130, 131, 135, 139, 149, 153, 157]

################################################################

def make_flat(): 

    util.mkdir(flat_dir)
    flat_num = np.arange(161,  177+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    calib.makeflat(scan_flat_frames, [], flat_dir + 'flat.fits', darks=False)

    return


def make_sky():

    util.mkdir(sky_dir)
    sky_num = np.arange(144, 148+1)
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
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame=flat_dir+"flat.fits")

    #EXPERIMENT 1:
    
    # Closed - LS_c
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c.fits'.format(ii) for ii in fnum_LS_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_LS_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - LS_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c.fits'.format(ii) for ii in fnum_LS_B2_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan.fits'.format(ii) for ii in fnum_LS_B2_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")
                
    # Closed - Mean_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c.fits'.format(ii) for ii in fnum_Mean_B2_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan.fits'.format(ii) for ii in fnum_Mean_B2_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    #EXPERIMENT 2:

    # Closed - nrej1
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej1_c.fits'.format(ii) for ii in fnum_nrej1]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej1_c_scan.fits'.format(ii) for ii in fnum_nrej1]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - nrej4
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej4_c.fits'.format(ii) for ii in fnum_nrej4]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej4_c_scan.fits'.format(ii) for ii in fnum_nrej4]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed - nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c.fits'.format(ii) for ii in fnum_nrej7]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_nrej7]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, \
                sky_frame=sky_dir + 'FLD2_sky.fits', flat_frame =flat_dir+"flat.fits")
    return



def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.find_stars(img_files, fwhm=10, threshold=20, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)

    # EXPERIMENT 1
                              
    # Closed - LS_c
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_LS_c]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)
    # Closed - LS_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan.fits'.format(ii) for ii in fnum_LS_B2_c]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)

    # Closed - Mean_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan.fits'.format(ii) for ii in fnum_Mean_B2_c]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)

    # EXPERIMENT 2

    # Closed - nrej1
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej1_c_scan.fits'.format(ii) for ii in fnum_nrej1]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)

    # Closed - nrej4
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej4_c_scan.fits'.format(ii) for ii in fnum_nrej4]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)

    # Closed - nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_nrej7]
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

    # EXPERIMENT 1:
    
    # Closed - LS_c
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_LS_c]
    stats_file = stats_dir + 'stats_closed_LS.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    # Closed - LS_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_B2_c_scan.fits'.format(ii) for ii in fnum_LS_B2_c]
    stats_file = stats_dir + 'stats_closed_LS_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - Mean_B2_c
    img_files = [data_dir + 'obj{0:03d}threeWFSMean_B2_c_scan.fits'.format(ii) for ii in fnum_Mean_B2_c]
    stats_file = stats_dir + 'stats_closed_Mean_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # EXPERIMENT 2

    # Closed - nrej1
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej1_c_scan.fits'.format(ii) for ii in fnum_nrej1]
    stats_file = stats_dir + 'stats_closed_nrej1.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - nrej4
    img_files = [data_dir + 'obj{0:03d}threeWFSLS_nrej4_c_scan.fits'.format(ii) for ii in fnum_nrej4]
    stats_file = stats_dir + 'stats_closed_nrej4.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)
    
    # Closed - nrej7
    img_files = [data_dir + 'obj{0:03d}threeWFS_LS_c_scan.fits'.format(ii) for ii in fnum_nrej7]
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
    
