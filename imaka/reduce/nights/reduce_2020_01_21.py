
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

#root_dir = '//Volumes/DATA5/imaka/20200121/sta/'
root_dir = '/g/lu/data/imaka/onaga/20200121/sta/'

obs_obj = 'Beehive-W' #observed object in the sky
sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + F'{obs_obj}/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + F'reduce/{obs_obj}/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

fnum_o = [3, 4, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, \
          32, 33, 34, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 59, 63, 65, 67, 68, 70, 72, 73, 74, 75, \
          76, 77, 79, 81, 83, 85, 87, 89, 91, 93, 107, 109, 111, 113, 115]

fnum_c =[8, 10, 36, 37, 39, 41, 43, 45, 47, 49, 51, 53, 57, 60, 61, 62, 64, 66, 69, 71, 78, 80, 82, 84, 86, 88, \
         90, 92, 106, 108, 110, 112, 114]

def make_flat(): 

    util.mkdir(flat_dir)
    
    flat_num = np.arange(221, 230+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    calib.makeflat(scan_flat_frames, [], flat_dir + 'flat.fits', darks=False)

    return


def make_sky():

#     util.mkdir(sky_dir)

#     sky_num = np.arange(0, 3)
#     sky_frames = ['{0:s}sta{1:03d}_o.fits'.format(root_dir+'HR2286/', ss) for ss in sky_num]
#     reduce_STA.treat_overscan(sky_frames)
#     scan_sky_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(root_dir+'HR2286/', ss) for ss in sky_num]
#     calib.makedark(scan_sky_frames, sky_dir+F'{obs_obj}_sky.fits')
    
    return


def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'sta{0:03d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'sta{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'beehive_sky_60.fits', flat_frame=None)

    # Closed - x10LS5WFS
    img_files = [data_dir + 'sta{0:03d}x10LS5WFS_c.fits'.format(ii) for ii in fnum_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'sta{0:03d}x10LS5WFS_c_scan.fits'.format(ii) for ii in fnum_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir + 'beehive_sky_60.fits', flat_frame =None)


    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    reduce_STA.find_stars(img_files, fwhm=8, threshold=10, min_flux = 10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=None, mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    
    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    reduce_STA.find_stars(img_files, fwhm=7, threshold=30, min_flux = 10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=None, mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=25)
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    stats_file = stats_dir + 'stats_open.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - x10LS5WFS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    stats_file = stats_dir + 'stats_closed_x10LS5WFS_c.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    return


def append_massdimm():

    massdimm.fetch_data('20200121', massdimm_dir)
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
    open_images = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'sta{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + F'{obs_obj}_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - x10LS5WFS
    closed_images = [out_dir + 'obj{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    closed_starlists = [out_dir + 'obj{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c]
    closed_output_root = stacks_dir + F'{obs_obj}_stack_x10LS5WFS_c'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    return


def analyze_stacks():

    open_img_files = [stacks_dir + F'{obs_obj}_stack_open.fits']

    closed_img_files = [stacks_dir + F'{obs_obj}_stack_x10LS5WFS_c.fits']
    
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
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.fourfilt(img_files, starlists)

    #Closed Loop - x10LS5WFS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.fourfilt(img_files, starlists)

    return

def calc_star_stats_fourfilt():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    

    # R
    stats_file = stats_dir + 'stats_open_R.fits'
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_R_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='R')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # V
    stats_file = stats_dir + 'stats_open_V.fits'
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_V_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='V')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # B
    stats_file = stats_dir + 'stats_open_B.fits'
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_B_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='B')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # I
    stats_file = stats_dir + 'stats_open_I.fits'
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_I_stars.txt'.format(ii) for ii in fnum_o]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='I')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    #Closed Loop - x10LS5WFS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]

    # R
    stats_file = stats_dir + 'stats_closed_R.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_R_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='R')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # V
    stats_file = stats_dir + 'stats_closed_V.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_V_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='V')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # B
    stats_file = stats_dir + 'stats_closed_B.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_B_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='B')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # I
    stats_file = stats_dir + 'stats_closed_I.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_I_stars.txt'.format(ii) for ii in fnum_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, filt='I')
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    return
