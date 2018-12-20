import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
#import scipy
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

root_dir = '//Volumes/DATA5/imaka/20181219/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Orion/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Orion/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

fnum_o = [36, 38, 41, 43, 45, 47, 49, 55, 58, 61, 64, 67, 73, 77, 80, 83, 86, 89, 96, 99, 101, 103, 106]
fnum_c_B2 = [35, 37, 40, 42, 44, 46, 48, 53, 56, 59, 62, 65, 68, 71, 74, 78, 81, 84, 87, 90, 94, 97]
fnum_c_n7 = [52, 54, 57, 60, 63, 66, 69, 70, 72, 75, 79, 82, 85, 88, 91, 93, 95, 98, 100, 105, 107]


def make_flat(): 

    util.mkdir(flat_dir)
    
    flat_num = np.arange(0,  9+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    calib.makeflat(scan_flat_frames, [], flat_dir + 'flat.fits', darks=False)

    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(14, 33+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'orion_sky.fits')
    
    return


def reduce_orion():

    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'obj{0:03d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o]
    areduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'orion_sky.fits', flat_frame=flat_dir+"flat.fits")

    # Closed Loop - B2
    img_files = [data_dir + 'obj{0:03d}LS_B2_c.fits'.format(ii) for ii in fnum_c_B2]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}LS_B2_c_scan.fits'.format(ii) for ii in fnum_c_B2]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'orion_sky.fits', flat_frame =flat_dir+"flat.fits")

    # Closed Loop - n7
    img_files = [data_dir + 'obj{0:03d}LSnrej7_Bin2_zc11_c.fits'.format(ii) for ii in fnum_c_n7]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'obj{0:03d}LS_Bin2_zc11_c_scan.fits'.format(ii) for ii in fnum_c_n7]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'orion_sky.fits', flat_frame =flat_dir+"flat.fits")

    return


def find_stars_orion():

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    #img_files = [out_dir + 'obj064_o_scan_clean.fits']
    reduce_fli.find_stars(img_files, fwhm=10, threshold=3, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
                              left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)
    
    #Closed Loop - LS
    img_files = [out_dir + 'obj{0:03d}LS_Bin2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    #img_files = [out_dir + 'obj062LS_c_scan_clean.fits']
    reduce_fli.find_stars(img_files, fwhm=7, threshold=3, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
                              left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)

    #Closed Loop - B2
    img_files = [out_dir + 'obj{0:03d}LS_c_scan_clean.fits'.format(ii) for ii in fnum_c_LS]
    #img_files = [out_dir + 'obj063LS_Bin2_c_scan_clean.fits']
    reduce_fli.find_stars(img_files, fwhm=7, threshold=3, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
                              left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    stats_file = stats_dir + 'stats_open.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - LS
    img_files = [out_dir + 'obj{0:03d}LS_c_scan_clean.fits'.format(ii) for ii in fnum_c_LS]
    stats_file = stats_dir + 'stats_closed_LS.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - B2
    img_files = [out_dir + 'obj{0:03d}LS_Bin2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    stats_file = stats_dir + 'stats_closed_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    return


def append_massdimm():

    massdimm.fetch_data('20181218', massdimm_dir)
    stats_tables = glob.glob(root_dir + 'reduce/stats/stats*.fits')

    for stats in stats_tables:
        if 'mdp.fits' not in stats:
            print('Adding MASS/DIMM to ' + stats)
            massdimm.append_mass_dimm(stats, massdimm_dir)
        else:
            print('Skipping ' + stats)

    return


def stack_orion():

    util.mkdir(stacks_dir)

    # Open Loop
    open_images = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'orion_stack_open'
    #reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - LS
    closed_images = [out_dir + 'obj{0:03d}LS_c_scan_clean.fits'.format(ii) for ii in fnum_c_LS]
    closed_starlists = [out_dir + 'obj{0:03d}LS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_LS]
    closed_output_root = stacks_dir + 'orion_stack_closed_LS'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - B2
    closed_images = [out_dir + 'obj{0:03d}LS_Bin2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    closed_starlists = [out_dir + 'obj{0:03d}LS_Bin2_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_B2]
    closed_output_root = stacks_dir + 'orion_stack_closed_BS'
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
    
