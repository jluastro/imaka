import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
#import scipy
import glob
from imaka.reduce import reduce_fli as redu
from imaka.reduce import calib
from imaka.reduce import util
from imaka.analysis import moffat
from astropy.stats import sigma_clipped_stats
import os, shutil
import pdb
from imaka.reduce import massdimm
from imaka.reduce import reduce_STA
import matplotlib
# matplotlib.use('Agg')

root_dir = '/g/lu/data/imaka/onaga/20210430/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Fld2/'
calib_dir = root_dir + 'reduce/calib/'
flat_dir = calib_dir # in case some name changes not caught
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
#old_cals_dir = root_dir + 'cals_from_20210403/'
massdimm_dir = root_dir + 'reduce/massdimm/'

# skies for bin2 are pulled from previous night
root_dir_2 = '/g/lu/data/imaka/onaga/20210429/sta/'
data_dir_bin2 = root_dir_2 + 'Fld2/'
sky_dir_bin2 = root_dir_2 + 'reduce/sky/'


# This night had two types of binning, suffix bin1, bin2 has been added to indicate

# Junk files -- see logs
# 

dict_suffix = {'open_bin1': '_o',
               'LS_bin1':   'LS_c',
               'docz_bin1': 'docz2_c',
               'open_bin2': '_o',
               'LS_bin2':   'LS_c',
               'docz_bin2': 'modal_c'}

dict_images = {'open_bin1': [68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94],
               'LS_bin1':   [67, 71, 75, 79, 83, 87, 91],
               'docz_bin1': [69, 73, 77, 81, 85, 89, 93],
               'open_bin2': [64, 66],
               'LS_bin2':   [63],
               'docz_bin2': [65]}

# These don't exist yet... make some temp files with zeros.
# have bin1 skies but not bin2, borrow from night2?
dict_skies = {'open_bin1':    'fld2_sky_180_bin1.fits',
              'LS_bin1':      'fld2_sky_180_bin1.fits',
              'docz_bin1':    'fld2_sky_180_bin1.fits',
              'open_bin2':    'fld2_sky_180_bin2.fits',
              'LS_bin2':      'fld2_sky_180_bin2.fits',
              'docz_bin2':    'fld2_sky_180_bin2.fits'}
    

def make_flat(): 
    """
    Make a flat... this will be with dome flats for now. 
    These are junk... replace with twilight flats. 
    """
    util.mkdir(flat_dir)

    # running for two different binnings
    flat_num_bin1 = np.arange(51, 54+1)
    flat_num_bin2 = np.arange(55, 62+1)
    flat_frames_bin1 = ['{0:s}twi{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num_bin1]
    flat_frames_bin2 = ['{0:s}twi{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num_bin2]
    reduce_STA.treat_overscan(flat_frames_bin1)
    reduce_STA.treat_overscan(flat_frames_bin2)
    
    scan_flat_frames_bin1 = ['{0:s}twi{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num_bin1]
    scan_flat_frames_bin2 = ['{0:s}twi{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num_bin2]
  
    calib.makeflat(scan_flat_frames_bin1, None, calib_dir + 'flat_bin1.fits', darks=False)
    calib.makeflat(scan_flat_frames_bin2, None, calib_dir + 'flat_bin2.fits', darks=False)

    # This mask tells us where not to search for stars.
    calib.make_mask(calib_dir + 'flat_bin1.fits', calib_dir + 'mask_bin1.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    calib.make_mask(calib_dir + 'flat_bin2.fits', calib_dir + 'mask_bin2.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=10, right_slice=10, top_slice=15, bottom_slice=15)

    return


def make_sky():

    util.mkdir(sky_dir)
    
    # bin1 skys were taken 0430
    sky_num_180_bin1 = np.arange(96, 103+1)
    sky_frames_180 = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num_180_bin1]
    reduce_STA.treat_overscan(sky_frames_180, remake=True)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num_180_bin1]
    calib.makedark(scan_sky_frames, sky_dir+'fld2_sky_180_bin1.fits')

    #bin2 skys were taken 0429, copy files
    shutil.copyfile(sky_dir_bin2+'fld2_sky_180.fits', sky_dir+'fld2_sky_180_bin2.fits')
    shutil.copyfile(sky_dir_bin2+'fld2_sky_180.list', sky_dir+'fld2_sky_180_bin2.list')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    # Loop through all the different data sets and reduce them.
    # for ke
    for key in dict_suffix.keys():
        print(key)
        bin = 'bin1' if 'bin1' in key else 'bin2'
	
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        reduce_fli.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=flat_dir+"flat_"+bin+".fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return


def find_stars_fld2():
    dict_fwhm = {'open_30': 7,
                 'LS_30':   3,
                 'docz_30': 3,
                 'open': 7,
                 'LS': 3,
                 'docz': 3,
                 'moda': 3}
    
    # Loop through all the different data sets and reduce them.
    # for key in dict_suffix.keys():
    for key in ['docz']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        reduce_fli.find_stars(img_files, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_flat=flat_dir+"flat_"+bin+".fits", mask_min=0.8, mask_max=1.4,
                              left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
                          
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
