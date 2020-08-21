
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

fnum_o = [3, 4, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 59, 63, 65, 67, 68, 70, 72, 73, 74, 75, \
          76, 77, 79, 81, 83, 85, 87, 89, 91, 93, 107, 109, 111, 113, 115]

fnum_c =[8, 36, 37, 39, 41, 43, 45, 47, 49, 51, 53, 57, 60, 61, 62, 64, 66, 69, 71, 78, 80, 82, 84, 86, 88, \
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

    util.mkdir(sky_dir)

    sky_num = np.arange(208, 220+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format('/g/lu/data/imaka/onaga/20200122/sta/Beehive-W/', ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format('/g/lu/data/imaka/onaga/20200122/sta/Beehive-W/', ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+F'{obs_obj}_sky.fits')
    
    return


def reduce_FLD2():
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'sta{0:03d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'sta{0:03d}_o_scan.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir+'Beehive-W_sky.fits', \
                            flat_frame=flat_dir+'flat.fits', fix_bad_pixels=True)

    # Closed - x10LS5WFS
    img_files = [data_dir + 'sta{0:03d}x10LS5WFS_c.fits'.format(ii) for ii in fnum_c]
    reduce_STA.treat_overscan(img_files)
    scan_img_files = [data_dir + 'sta{0:03d}x10LS5WFS_c_scan.fits'.format(ii) for ii in fnum_c]
    reduce_fli.clean_images(scan_img_files, out_dir, rebin=1, sky_frame=sky_dir+'Beehive-W_sky.fits', \
                            flat_frame=flat_dir+'flat.fits', fix_bad_pixels=True)


    return


def find_stars_FLD2():
    flat = flat_dir+'flat.fits'
    # Open Loop
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    reduce_STA.find_stars(img_files, fwhm=8, threshold=5, min_flux =13, min_x_fwhm=2.5, min_y_fwhm=2.5, N_passes=2, \
                          plot_psf_compare=False, mask_flat=None)
#                               mask_flat=flat, mask_min=0.752, mask_max=1.4, \
#                               left_slice =1247, right_slice=3518, top_slice=2646, bottom_slice=0)
    
    #Closed Loop - threeWFS_LS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    reduce_STA.find_stars(img_files, fwhm=8, threshold=10, min_flux =13, min_x_fwhm=2.0, min_y_fwhm=2.0, N_passes=2, \
                          plot_psf_compare=False, mask_flat=None)
#                               mask_flat=flat, mask_min=0.752, mask_max=1.4, \
#                               left_slice =1247, right_slice=3518, top_slice=2646, bottom_slice=0)
                          
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
    reduce_STA.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - x10LS5WFS
    closed_images = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c]
    closed_starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c]
    closed_output_root = stacks_dir + F'{obs_obj}_stack_x10LS5WFS_c'
    reduce_STA.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    return


def analyze_stacks():

    open_img_files = [stacks_dir + F'{obs_obj}_stack_open.fits']

    closed_img_files = [stacks_dir + F'{obs_obj}_stack_x10LS5WFS_c.fits']
    
    #Find stars in image
    reduce_STA.find_stars(open_img_files, fwhm=10, threshold=100, N_passes=2, plot_psf_compare=False, \
                          mask_flat=None)
#                               mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
#                               left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    reduce_STA.find_stars(closed_img_files, fwhm=7, threshold=30, N_passes=2, plot_psf_compare=False, \
                          mask_flat=None)
#                               mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
#                               left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
        
    # Calc stats on all the stacked images
    #reduce_fli.calc_star_stats(open_img_files+closed_img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return

###############################################
############ FOUR FILTER REDUCTION ############
###############################################

"""
Notes on rotation from log:
POS 1
at frame 61:  I(NW), V(SW), B(SE), R(NE)  (IRBV)
POS 2
at frame 82:  V(NW), I(NE), R(SE), B(SW)  (VIRB)
POS 3
at frame 106: B(NW), V(NE), I (SE), R(SW) (BVIR)
"""

# Open
rot_1_o = [65, 67, 70, 79, 81]
rot_2_o = [83, 85, 87, 89, 91, 93]
rot_3_o = [107, 109, 111, 113, 115]

# Closed
rot_1_c = [61, 64, 66, 69, 78, 80] 
rot_2_c = [82, 84, 86, 88, 90, 92] 
rot_3_c = [106, 108, 110, 112, 114]

rot_o_4 = rot_1_o + rot_2_o + rot_3_o
rot_c_4 = rot_1_c + rot_2_c + rot_3_c
    
def split_filt():

    # Rotation 1: open
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_1_o]
    reduce_STA.four_filt_split(starlists, 'IRBV')
    
    # Rotation 2: open
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_2_o]
    reduce_STA.four_filt_split(starlists, 'VIRB')

    # Rotation 3: open
    starlists = [out_dir + 'sta{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_3_o]
    reduce_STA.four_filt_split(starlists, 'BVIR')

    # Rotation 1: closed
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_1_c]
    reduce_STA.four_filt_split(starlists, 'IRBV')
    
    # Rotation 2: closed
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_2_c]
    reduce_STA.four_filt_split(starlists, 'VIRB')

    # Rotation 3: closed
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.four_filt_split(starlists, 'BVIR')

    return

def calc_star_stats_fourfilt():
    util.mkdir(stats_dir)

    # Open Loop
    img_files = [out_dir + 'sta{0:03d}_o_scan_clean.fits'.format(ii) for ii in rot_o_4]
    

    # # Open Loop - B
    starlists =  [out_dir + 'sta{0:03d}_o_scan_clean_B_IRBV_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_B_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_B_BVIR_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_B.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

#     # Open Loop - V
    starlists =  [out_dir + 'sta{0:03d}_o_scan_clean_V_IRBV_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_V_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_V_BVIR_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_V.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

#     # Open Loop - R
    starlists =  [out_dir + 'sta{0:03d}_o_scan_clean_R_IRBV_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_R_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_R_BVIR_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_R.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

#     # Open Loop - I
    starlists =  [out_dir + 'sta{0:03d}_o_scan_clean_I_IRBV_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_I_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'sta{0:03d}_o_scan_clean_I_BVIR_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_I.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    #Closed Loop - x10LS5WFS
    img_files = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean.fits'.format(ii) for ii in rot_c_4]

    # R
    stats_file = stats_dir + 'stats_closed_R.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_R_IRBV_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_R_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_R_BVIR_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # V
    stats_file = stats_dir + 'stats_closed_V.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_V_IRBV_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_V_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_V_BVIR_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # B
    stats_file = stats_dir + 'stats_closed_B.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_B_IRBV_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_B_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_B_BVIR_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # I
    stats_file = stats_dir + 'stats_closed_I.fits'
    starlists = [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_I_IRBV_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_I_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'sta{0:03d}x10LS5WFS_c_scan_clean_I_BVIR_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    return
