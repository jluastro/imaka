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

root_dir = '/g/lu/data/imaka/onaga/20200122/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Beehive-W/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Beehive-W/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

fnum_o = [113, 118, 125, 130, 137, 142, 147, 152, 157, 162, 165, 167, 172, 177, 182, 187, 192, 197]       # open loop
fnum_o_150 = [199, 202, 205]       # open loop, long exposure
fnum_c_5W = [132, 136, 146, 156, 168, 173, 178, 183, 190, 195] # cmat 0 (dropped 114 and 121.. failed offsets load and bad filename)
fnum_c_5W_150 = [198, 200, 201, 203, 204]       # cmat 0, long exposure
fnum_c_4W = [115, 122, 135, 145, 155] # cmat 1 -- first two: slope offsets failed to load.
fnum_c_3W = [116, 123, 134, 144, 154] # cmat 2
fnum_c_1W = [117, 124, 133, 143, 153] # cmat 3
fnum_c_cri2 = [119, 126, 141, 151, 161]  # cmat 4
fnum_c_cri3 = [120, 127, 140, 150, 160, 166, 170, 175, 181, 186, 191, 196]      # cmat 5
fnum_c_cri4 = [163, 171, 176, 180, 185, 189, 194]      # cmat 5
fnum_c_9zer5 = [128, 139, 149, 159, 164, 169, 174, 179, 184, 188, 193]     # cmat 6
fnum_c_9zer4 = [129, 138, 148, 158]     # cmat 7

dict_suffix = {'open': '_o',
               '5W': 'LS5WFS_c',
               '4W': 'LS4WFS_c',
               '3W': 'LS3WFS_c',
               '1W': 'LS1WFS_c',
               'cri2': 'docrime2_c',
               'cri3': 'docrime3_c',
               'cri4': 'docrime4_c',
               '9zer5': '_9zer5WFS_c',
               '9zer4': '_9zerWFS4_c',
               'open_long': '_o',                       
               '5W_long': 'LS5WFS_c'}
    
dict_images = {'open': fnum_o,
               '5W': fnum_c_5W,
               '4W': fnum_c_4W,
               '3W': fnum_c_3W,
               '1W': fnum_c_1W,
               'cri2': fnum_c_cri2,
               'cri3': fnum_c_cri3,
               'cri4': fnum_c_cri4,
               '9zer5': fnum_c_9zer5,
               '9zer4': fnum_c_9zer4,
               'open_long': fnum_o_150,
               '5W_long': fnum_c_5W_150}

dict_skies =  {'open': 'beehive_sky_60.fits',
               '5W': 'beehive_sky_60.fits',
               '4W': 'beehive_sky_60.fits',
               '3W': 'beehive_sky_60.fits',
               '1W': 'beehive_sky_60.fits',
               'cri2': 'beehive_sky_60.fits',
               'cri3': 'beehive_sky_60.fits',
               'cri4': 'beehive_sky_60.fits',
               '9zer5': 'beehive_sky_60.fits',
               '9zer4': 'beehive_sky_60.fits',
               'open_long': 'beehive_sky_150.fits',
               '5W_long': 'beehive_sky_150.fits'}
    

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

    sky_num = np.arange(208, 217+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'beehive_sky_60.fits')

    sky_num = np.arange(218, 220+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'beehive_sky_150.fits')
    
    return


def reduce_beehive():

    util.mkdir(out_dir)


    # Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
                  
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        reduce_fli.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=flat_dir + "flat.fits")

    return


def find_stars_beehive():

    # Open Loop
    img_files = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    #reduce_fli.find_stars(img_files, fwhm=8, threshold=6, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
    #                          left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)
    
    #Closed Loop - 4W
    img_files = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c_4W]
    #reduce_fli.find_stars(img_files, fwhm=7, threshold=6, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
    #                          left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)

    #Closed Loop - B2
    img_files = [out_dir + 'obj{0:03d}LS4WFS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    #reduce_fli.find_stars(img_files, fwhm=7, threshold=6, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
    #                          left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)

    #Closed Loop - zc
    img_files = [out_dir + 'obj{0:03d}LS4WFS_zc11_c_scan_clean.fits'.format(ii) for ii in fnum_c_zc]
    #reduce_fli.find_stars(img_files, fwhm=7, threshold=6, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.8, mask_max=1.4, \
    #                          left_slice =20, right_slice=20, top_slice=25, bottom_slice=25)

    #Tip tilt
    img_files = [out_dir + 'obj{0:03d}tip_tilt_scan_clean.fits'.format(ii) for ii in fnum_tt]
    reduce_fli.find_stars(img_files, fwhm=7, threshold=6, N_passes=2, plot_psf_compare=False, \
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

    #Closed Loop - B2
    img_files = [out_dir + 'obj{0:03d}LS4WFS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    stats_file = stats_dir + 'stats_closed_B2.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - 4W
    img_files = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c_4W]
    stats_file = stats_dir + 'stats_closed_4W.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #Closed Loop - n1
    img_files = [out_dir + 'obj{0:03d}LS4WFS_zc11_c_scan_clean.fits'.format(ii) for ii in fnum_c_zc]
    stats_file = stats_dir + 'stats_closed_zc.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)

    #tip tilt                                                          
    img_files = [out_dir + 'obj{0:03d}tip_tilt_scan_clean.fits'.format(ii)\
 for ii in fnum_tt]
    stats_file = stats_dir + 'stats_closed_tt.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
    moffat.fit_moffat(img_files, stats_file)


    return


def append_massdimm():

    massdimm.fetch_data('20181223', massdimm_dir)
    stats_tables = glob.glob(root_dir + 'reduce/stats/stats*.fits')

    for stats in stats_tables:
        if 'mdp.fits' not in stats:
            print('Adding MASS/DIMM to ' + stats)
            massdimm.append_mass_dimm(stats, massdimm_dir)
        else:
            print('Skipping ' + stats)

    return


def stack_beehive():

    util.mkdir(stacks_dir)

    # Open Loop
    open_images = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'beehive_stack_open'
    #reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - 4W
    closed_images = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c_4W]
    closed_starlists = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_4W]
    closed_output_root = stacks_dir + 'beehive_stack_closed_4W'
    #reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed Loop - B2
    closed_images = [out_dir + 'obj{0:03d}LS4WFS_B2_c_scan_clean.fits'.format(ii) for ii in fnum_c_B2]
    closed_starlists = [out_dir + 'obj{0:03d}LS4WFS_B2_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_B2]
    closed_output_root = stacks_dir + 'beehive_stack_closed_B2'
    #reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    #Closed Loop - zc
    closed_images = [out_dir + 'obj{0:03d}LS4WFS_zc11_c_scan_clean.fits'.format(ii) for ii in fnum_c_zc]
    closed_starlists = [out_dir + 'obj{0:03d}LS4WFS_zc11_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_zc]
    closed_output_root = stacks_dir + 'beehive_stack_closed_zc'
    #reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    #tt
    closed_images = [out_dir + 'obj{0:03d}tip_tilt_scan_clean.fits'.format(ii) for ii in fnum_tt]
    closed_starlists = [out_dir + 'obj{0:03d}tip_tilt_scan_clean_stars.txt'.format(ii) for ii in fnum_tt]
    closed_output_root = stacks_dir + 'beehive_stack_tiptilt'
    #reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return


def analyze_stacks():

    open_img_files = [stacks_dir + 'beehive_stack_open.fits']

    closed_img_files = [#stacks_dir + 'beehive_stack_closed_4W.fits', \
                        stacks_dir + 'beehive_stack_closed_B2.fits', \
                        stacks_dir + 'beehive_stack_closed_zc.fits', \
                        stacks_dir + 'beehive_stack_tiptilt.fits']
    
    #Find stars in image
    #reduce_fli.find_stars(open_img_files, fwhm=10, threshold=10, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
    #                          left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    reduce_fli.find_stars(closed_img_files, fwhm=7, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
        
    # Calc stats on all the stacked images
    #reduce_fli.calc_star_stats(open_img_files+closed_img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return

###############################################
############ FOUR FILTER REDUCTION ############
###############################################

"""
Notes on rotation from log:
at frame 16: B (NW), V(NE), R(SW), I(SE)
at frame 49: V (NW), B(SW), I(NE), R(SE)
at frame 75: B (NW), V(NE), R(SW), I(SE)
at frame 81: B (SE), V(SW), R(NE), I(NW)
at frame 99: switched filter to I band 
"""
    
rot_1_o = [30, 31, 32, 34, 36, 38, 40, 42, 44, 46, 48]
rot_1_o += [77, 80]
rot_2_o = [52, 53, 57, 59, 62, 65, 68, 71, 74]
rot_3_o = [83, 86, 89, 92, 95, 98]

rot_1_c = [33, 35, 37, 39, 41, 43, 45, 47] 
rot_1_c += [76, 79]
rot_2_c = [55, 58, 61, 64, 67, 70, 73] 
rot_3_c = [82, 85, 88, 91, 94, 97]

rot_o_4 = rot_1_o + rot_2_o + rot_3_o
rot_c_4 = rot_1_c + rot_2_c + rot_3_c


def split_filters():
   
    # Rotation 1: open
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_1_o]
    reduce_STA.four_filt_split(starlists, 'BVIR')
    
    # Rotation 2: open
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_2_o]
    reduce_STA.four_filt_split(starlists, 'VIRB')

    # Rotation 3: open
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in rot_3_o]
    reduce_STA.four_filt_split(starlists, 'IRBV')

    # Rotation 1: closed
    starlists = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_1_c]
    reduce_STA.four_filt_split(starlists, 'BVIR')
    
    # Rotation 2: closed
    starlists = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_2_c]
    reduce_STA.four_filt_split(starlists, 'VIRB')

    # Rotation 3: closed
    starlists = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_stars.txt'.format(ii) for ii in rot_3_c]
    reduce_STA.four_filt_split(starlists, 'IRBV')
    
    return
    

def calc_fourfilt_stats():
    
    # Open Loop - B
    img_files =  [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in rot_o_4]
    starlists =  [out_dir + 'obj{0:03d}_o_scan_clean_B_BVIR_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_B_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_B_IRBV_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_B.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # Open Loop - V
    img_files =  [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in rot_o_4]
    starlists =  [out_dir + 'obj{0:03d}_o_scan_clean_V_BVIR_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_V_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_V_IRBV_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_V.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # Open Loop - R
    img_files =  [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in rot_o_4]
    starlists =  [out_dir + 'obj{0:03d}_o_scan_clean_R_BVIR_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_R_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_R_IRBV_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_R.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # Open Loop - I
    img_files =  [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in rot_o_4]
    starlists =  [out_dir + 'obj{0:03d}_o_scan_clean_I_BVIR_stars.txt'.format(ii) for ii in rot_1_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_I_VIRB_stars.txt'.format(ii) for ii in rot_2_o]
    starlists += [out_dir + 'obj{0:03d}_o_scan_clean_I_IRBV_stars.txt'.format(ii) for ii in rot_3_o]
    stats_file = stats_dir + 'stats_open_I.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # Closed Loop - B
    img_files =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in rot_c_4]
    starlists =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_B_BVIR_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_B_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_B_IRBV_stars.txt'.format(ii) for ii in rot_3_c]
    stats_file = stats_dir + 'stats_closed_B.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # Closed Loop - V
    img_files =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in rot_c_4]
    starlists =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_V_BVIR_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_V_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_V_IRBV_stars.txt'.format(ii) for ii in rot_3_c]
    stats_file = stats_dir + 'stats_closed_V.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    # Closed Loop - R
    img_files =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in rot_c_4]
    starlists =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_R_BVIR_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_R_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_R_IRBV_stars.txt'.format(ii) for ii in rot_3_c]
    stats_file = stats_dir + 'stats_closed_R.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)
    
    # Closed Loop - I
    img_files =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in rot_c_4]
    starlists =  [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_I_BVIR_stars.txt'.format(ii) for ii in rot_1_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_I_VIRB_stars.txt'.format(ii) for ii in rot_2_c]
    starlists += [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_I_IRBV_stars.txt'.format(ii) for ii in rot_3_c]
    stats_file = stats_dir + 'stats_closed_I.fits'
    reduce_STA.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
    moffat.fit_moffat(img_files, stats_file, starlists=starlists)

    return

def stack_beehive_I():

    # I Band - open
    images = [out_dir + 'obj{0:03d}_o_scan_clean.fits'.format(ii) for ii in fnum_o_I]
    starlists = [out_dir + 'obj{0:03d}_o_scan_clean_stars.txt'.format(ii) for ii in fnum_o_I]
    output_root = stacks_dir + 'beehive_stack_open_I'
    reduce_fli.shift_and_add(images, starlists, output_root, method='mean')

    # I Band - closed
    images = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean.fits'.format(ii) for ii in fnum_c_I]
    starlists = [out_dir + 'obj{0:03d}LS4WFS_c_scan_clean_stars.txt'.format(ii) for ii in fnum_c_I]
    output_root = stacks_dir + 'beehive_stack_closed_I'
    reduce_fli.shift_and_add(images, starlists, output_root, method='mean')

    return

def analyze_stacks_I():

    open_img_files = [stacks_dir + 'beehive_stack_open_I.fits']

    closed_img_files = [stacks_dir + 'beehive_stack_closed_I.fits']
    
    #Find stars in image
    #reduce_fli.find_stars(open_img_files, fwhm=10, threshold=10, N_passes=2, plot_psf_compare=False, \
    #                          mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
    #                          left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    reduce_fli.find_stars(closed_img_files, fwhm=7, threshold=10, N_passes=2, plot_psf_compare=False, \
                              mask_flat=flat_dir+"flat.fits", mask_min=0.7, mask_max=1.4, \
                              left_slice =25, right_slice=0, top_slice=25, bottom_slice=0)
    
    return
