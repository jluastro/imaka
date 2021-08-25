## reduce_2021_07_24.py
##########################
## edited by Eden McEwen
## July 2021
## NOTE: no science images for this night
## Primarily flats and darks testing

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

root_dir = '/g/lu/data/imaka/onaga/20210724/sta/'

data_dir = root_dir + 'Fld2/'
twi_dir = root_dir + 'twilights/'
cal_dir = root_dir + 'calunit/'
dome_dir = root_dir + 'dome/'
dark_dir = root_dir + 'dark/'

sky_dir = root_dir + 'reduce/sky/' 
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
massdimm_dir = root_dir + 'reduce/massdimm/'

## Junk files -- see logs

dict_suffix = {'open': '_o',
               'LS': 'LS_c',
               'docz': 'docz2_c',
               'doczskycl': 'doczskycl_c',
               'z10glao': 'z10glao_c',
               'z10scao': 'Z10scao_c'}

dict_images = {'open':      [],
               'LS':        [],
               'docz':      [],
               'doczskycl': [],
               'z10glao':   [],
               'z10scao':   []
              }

dict_skies = {'open':      'fld2_sky.fits',
              'LS':        'fld2_sky.fits',
              'docz':      'fld2_sky.fits',
              'doczskycl': 'fld2_sky.fits',
              'z10glao': 'fld2_sky.fits',
              'z10scao': 'fld2_sky.fits'
              }

dict_fwhm = {'open': 12,
             'LS': 5,
             'docz': 5,
             'doczskycl': 5,
             'z10glao': 5,
             'z10scao': 5
            }  

###############################################
### REDUCTION
###############################################


def make_flat(): 
    """
    Makes flat and data mask. 
    Just for single filter I band
    """
    util.mkdir(calib_dir)
    
    ## Copy flat from a previous night
    shutil.copyfile(root_dir + '../../20210723/sta/reduce/calib/domeflat_I.fits', calib_dir + 'flat_I.fits')
    
    ## Creating flat from range, I band only
    #flat_num = np.arange(37, 49+1)
    #flat_frames = ['{0:s}dome_{1:03d}.fits'.format(dome_dir, ss) for ss in flat_num]
    #scan_flat_frames = ['{0:s}dome_{1:03d}_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    
    #reduce_STA.treat_overscan(flat_frames)
    #calib.makeflat(scan_flat_frames, None, calib_dir + 'domeflat_I.fits', darks=False)

    ## Make a mask to use when calling find_stars.
    calib.make_mask(calib_dir + 'domeflat_I.fits', calib_dir + 'mask.fits',
                       mask_min=0.5, mask_max=1.8,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(97, 101+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    
    reduce_STA.treat_overscan(sky_frames)
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky.fits')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    #for key in ['set_name']: ## Single key setup
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
        redu.clean_images(scn_files, out_dir, rebin=1,
                                    sky_frame=sky_dir + sky,
                                    flat_frame=calib_dir + "flat.fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return

###############################################
### ANALYSIS
###############################################

def find_stars_fld2():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]

        redu.find_stars(img_files, fwhm=fwhm, threshold=6, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask.fits')
        
    ## DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)
    
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Catalog: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        stats_file = stats_dir + 'stats_' + key + '.fits'
        
        redu.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0])
    # stats_file = stats_dir + 'stats_LS_c.fits'
    # redu.calc_star_stats(image_file, stats_file, flux_percent=0.2)
        
    return


def append_massdimm():
    date_str_start = root_dir.index('20')
    date_str = root_dir[date_str_start:date_str_start+8]
    print(f'Fetching MASS/DIMM for {date_str}')

    massdimm.fetch_data(date_str, massdimm_dir)
    stats_tables = glob.glob(root_dir + 'reduce/stats/stats*.fits')

    for stats in stats_tables:
        if 'mdp.fits' not in stats:
            print('Adding MASS/DIMM to ' + stats)
            try:
                massdimm.append_mass_dimm(stats, massdimm_dir)
            except Exception as e:
                print(" => ERROR with " + stats)
                print(" => ", e)
        else:
            print('Skipping ' + stats)

    return


def stack():

    util.mkdir(stacks_dir)

    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        output_root = stacks_dir + 'fld2_stack_' + suf
        
        redu.shift_and_add(img_files, starlists, output_root, method='mean')
        
    return


def analyze_stacks():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    all_images = []
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False, mask_file=calib_dir + 'mask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return


###############################################
############ FOUR FILTER REDUCTION ############
###############################################

"""
Notes on rotation from log:
POS 1
at frame 61:  R(NW), V(SW), B(SE), R(NE)  (IRBV)
POS 2
at frame 82:  V(NW), I(NE), R(SE), B(SW)  (VIRB)
POS 3
at frame 106: B(NW), V(NE), I (SE), R(SW) (BVIR)
"""

## Splitting the initial dictonary keys above means they don't need to be again split here.
## Intentionally redundant here

filters = ['B', 'V', 'R', 'I']

# Open
rot_1_o = []      # key: o_4F_1
rot_2_o = []  # key: o_4F_2
rot_3_o = [] # key: o_4F_3

# Closed
rot_1_c = []  # key: c_4F_1
rot_2_c = []  # key: c_4F_2
rot_3_c = [] # key: c_4F_3

rot_o_4 = rot_1_o + rot_2_o + rot_3_o
rot_c_4 = rot_1_c + rot_2_c + rot_3_c

# files
dict_images_rot = {'o_4F_1': rot_1_o,
                    'o_4F_2': rot_2_o,
                    'o_4F_3': rot_3_o,
                    'c_4F_1': rot_1_c,
                    'c_4F_2': rot_2_c,
                    'c_4F_3': rot_3_c}

# suffixes
dict_suffix_rot = {'o_4F_1': '_o',
                  'o_4F_2': '_o',
                  'o_4F_3': '_o',
                  'c_4F_1': 'x10LS5WFS_c',
                  'c_4F_2': 'x10LS5WFS_c',
                  'c_4F_3': 'x10LS5WFS_c'}

# Filter Order
dict_orders_rot = {'o_4F_1': 'IRBV',
                  'o_4F_2': 'VIRB',
                  'o_4F_3': 'BVIR',
                  'c_4F_1': 'IRBV',
                  'c_4F_2': 'VIRB',
                  'c_4F_3': 'BVIR'}

def make_flat_filter(): 
    """
    Makes flat and data mask. 
    For four filter, dome flat was made with two integration times, 20(2) for BVR 60(45) for I
    """
    util.mkdir(calib_dir)
    
    ## I quad flat (60)
    dark_num = np.arange(108, 113+1)
    flat_num = np.arange(130, 135+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    flat_frames = ['{0:s}dome_{1:03d}.fits'.format(dome_dir, ss) for ss in flat_num]
    scan_flat_frames = ['{0:s}dome_{1:03d}_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    
    reduce_STA.treat_overscan(dark_frames)
    reduce_STA.treat_overscan(flat_frames)
    #calib.makeflat(scan_flat_frames, scan_dark_frames, calib_dir + 'domeflat_60.fits', darks=True, fourfilter=True)
    
    ## BVR quad flat (20)
    dark_num = np.arange(118, 123+1)
    flat_num = np.arange(124, 129+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    flat_frames = ['{0:s}dome_{1:03d}.fits'.format(dome_dir, ss) for ss in flat_num]
    scan_flat_frames = ['{0:s}dome_{1:03d}_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    
    reduce_STA.treat_overscan(dark_frames)
    reduce_STA.treat_overscan(flat_frames)
    #calib.makeflat(scan_flat_frames, scan_dark_frames, calib_dir + 'domeflat_20.fits', darks=True, fourfilter=True)
    
    # Combining two flats based on filter orientation
    
    
    
    ## Make a mask to use when calling find_stars.
    #calib.make_mask(calib_dir + 'domeflat_60.fits', calib_dir + 'mask_60.fits',
    #                   mask_min=0.5, mask_max=1.8,
    #                   left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    return


def split_filters():
    # Split all starlists by filter, given rotation order
    for key in dict_suffix_rot.keys():
        img = dict_images_rot[key]
        suf = dict_suffix_rot[key]
        odr = dict_orders_rot[key]
        
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        reduce_STA.four_filt_split(starlists, odr)
    
    return
    

def calc_fourfilt_stats():   
    # Getting unique suffixes (loop states):
    suffixes = list(set(dict_suffix_rot.values()))
    
    # Grouping by suffixes 
    for suf in suffixes:
        # keys with given suffix
        keys = [key for key in dict_suffix_rot.keys() if dict_suffix_rot[key] == suf]
        
        # Iterating through filters
        for f in filters:
            stats_file = stats_dir + f'stats_{suf}_{f}.fits'
            img_files = []
            starlists = []
            
            for key in keys:
                img = dict_images_rot[key]
                odr = dict_orders_rot[key]
                
                img_files += [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
                starlists += [out_dir + 'sta{img:03d}{suf:s}_scan_clean_{f:s}_{odr:s}_stars.txt'.format(img=ii, suf=suf, f=f, odr=odr) for ii in img]
            
            print(f"Calc Star_Stats: {suf} Filter: {f}")
            reduce_fli.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
            print("Starting moffat fitting")
            #moffat.fit_moffat(img_files, stats_file, starlists=starlists)
            
            ## DEBBUG: SINGLE THREAD
            # reduce_fli.calc_star_stats_single(img_files[0], starlists[0], True)
    
    return



