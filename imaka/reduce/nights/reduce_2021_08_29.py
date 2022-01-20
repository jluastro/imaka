## reduce_2021_08_29.py
##########################
## edited by Eden McEwen
## August 2021
## A four filter run
## Did I Band, donuts,  position 3 (VBRI), position 1 (RIVB) this night

import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
#import scipy
import glob
from imaka.reduce import reduce_fli as redu
from imaka.reduce import reduce_STA
from imaka.reduce import calib
from imaka.reduce import util
from imaka.analysis import moffat
from astropy.stats import sigma_clipped_stats
import os, shutil
import pdb
from imaka.reduce import massdimm
import matplotlib
# matplotlib.use('Agg')

night = '20210829'
root_dir = f'/g/lu/data/imaka/onaga/{night}/sta/'

data_dir = root_dir + 'Fld2/'
twi_dir = root_dir + 'twilights/'
cal_dir = root_dir + 'calunit/'
dome_dir = root_dir + 'dome/'
dark_dir = root_dir + 'dark/'
screen_dir = root_dir + 'Screen/'

sky_dir = root_dir + 'reduce/sky/' 
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
massdimm_dir = root_dir + 'reduce/massdimm/'

## pulling darks from a different date
alt_root_dir =f'/g/lu/data/imaka/onaga/20210724/sta/'
alt_dark_dir = alt_root_dir + 'dark/'

## Junk files -- see logs
## 

dict_suffix = {'LS_I':   'LS_c',
               'docz_I': 'docz2_c',
               'open_I': '_o',
               'tt_I':   'tt_c', 
               'LS_VBRI':   'LS_c',
               'docz_VBRI': 'docz2_c',
               'open_VBRI': '_o',
               'tt_VBRI':   'tt_c',
               'LS_RIVB':   'LS_c',
               'docz_RIVB': 'docz2_c',
               'open_RIVB': '_o',
               'tt_RIVB':   'tt_c'
              }

dict_images = {'LS_I':    [6, 7, 11, 15, 19, 23],
               'docz_I':  [8, 12, 16, 20, 24],
               'open_I':  [9, 13, 17, 21, 25],
               'tt_I':    [10, 14, 18, 22, 26],
               'LS_VBRI':    [35, 39, 40, 44],
               'docz_VBRI':  [36, 41, 45],
               'open_VBRI':  [37, 42, 46],
               'tt_VBRI':    [38, 43, 47],
               'LS_RIVB':   [55, 59, 63, 67],
               'docz_RIVB': [56, 60, 64, 68],
               'open_RIVB': [57, 61, 65, 69],
               'tt_RIVB':   [58, 62, 66, 70]
              }

dict_fwhm = {'LS_I': 6,
             'docz_I': 6,
             'open_I': 12,
             'tt_I': 6,
             'LS_VBRI': 6,
             'docz_VBRI': 6,
             'open_VBRI': 12,
             'tt_VBRI':  6,
             'LS_VBRI': 6,
             'docz_VBRI': 6,
             'open_VBRI': 12,
             'tt_VBRI':  6
            }

dict_filt = {'LS_I':  'I',
            'docz_I': 'I',
            'open_I': 'I',
            'tt_I':   'I', 
             'LS_VBRI':   'VBRI',
             'docz_VBRI': 'VBRI',
             'open_VBRI': 'VBRI',
             'tt_VBRI':   'VBRI',
             'LS_RIVB':   'RIVB',
             'docz_RIVB': 'RIVB',
             'open_RIVB': 'RIVB',
             'tt_RIVB':   'RIVB'
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
    
    ## Copy flat taken the a previous night
    # shutil.copyfile(root_dir + '../../20210828/sta/reduce/calib/flat_1p_VBRI.fits', calib_dir + 'flat_VBRI.fits')
    
    ## Creating flat from range, I band only
    flat_num = np.arange(1, 5+1)
    flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(screen_dir, ss) for ss in flat_num]
    scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(screen_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    calib.makeflat(scan_flat_frames, None, calib_dir + 'flat_I.fits', darks=False)

    ## Make a mask to use when calling find_stars.
    calib.make_mask(calib_dir + 'flat_I.fits', calib_dir + 'mask_I.fits',
                       mask_min=0.5, mask_max=1.8,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    #calib.make_mask(calib_dir + 'flat_VBRI.fits', calib_dir + 'mask_VBRI.fits',
    #                   mask_min=0.5, mask_max=1.8,
    #                   left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    return


def make_sky():

    util.mkdir(sky_dir)
    
    ## I Band
    ## NO SKY TAKEN => use a 180 dark
    dark_num = np.arange(103, 111+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames =  ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    reduce_STA.treat_overscan(dark_frames)
    calib.makedark(scan_dark_frames, sky_dir + 'fld2_sky_I.fits')
    
    ## VBRI pos 3
    sky_num = np.arange(48, 54)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky_VBRI.fits')
    
    ## RIVB pos 1
    sky_num = np.arange(71, 77+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky_RIVB.fits')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    #for key in ['tt_IVBR', 'LS_VBRI', 'open_VBRI', 'docz_VBRI', 'tt_VBRI']: ## Single key setup
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]
        filt = dict_filt[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('   Filter: ', filt)
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        redu.clean_images(scn_files, out_dir, rebin=1,
                                    sky_frame=sky_dir + f"fld2_sky_{filt}.fits",
                                    flat_frame=calib_dir + f"flat_{filt}.fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return

###############################################
### ANALYSIS
###############################################

def find_stars_fld2():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    #for key in ['set_name']:
    #for key in dict_suffix.keys():
        
    #    img = dict_images[key]
    #    suf = dict_suffix[key]
    #    filt = dict_filt[key]
    #    fwhm = dict_fwhm[key]

    #    print('Working on: {1:s}  {0:s}'.format(key, suf))
    #    print('   Images: ', img)
    #    print('      filt: ', filt)
        
    #    img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]

     #   redu.find_stars(img_files, fwhm=fwhm, threshold=8, N_passes=2, plot_psf_compare=False, mask_file=calib_dir+f'mask_{filt}.fits')
        
    ## DEBUG - single threaded
    key_test = 'LS_I'
    fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    image_file = fmt.format(dir=out_dir, img=dict_images[key_test][0], suf=dict_suffix[key_test]) 
    redu.find_stars_single(image_file, dict_fwhm[key_test], 3, 2, False, calib_dir + 'mask_I.fits')
                          
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
        filt = dict_filt[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False, mask_file=calib_dir + f'mask_{filt}.fits')

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
Notes on rotation from 07_23 and 07_24 log:
POS 1: RIVB - R(NW), I(NE), V(SE), B(SW). 
POS 2: IVBR - I(NW), V(NE), B(SE), R(SW). 
POS 3: VBRI - V(NW), B(NE), R(SE), I(SW). 
POS 4: BRIV - B(NW), R(NE), I(SE), V(SW).  

"""

## Splitting the initial dictonary keys above means they don't need to be again split here.

filters = [ 'V', 'R', 'I', 'B']

## Catching any missed name changes
dict_orders_rot = dict_filt
dict_images_rot = dict_images
dict_suffix_rot = dict_suffix


def make_flat_filter(): 
    """
    Makes flat and data mask. 
     - makeflat(fourfilter=True) normalizes each quadrant independently
     - combine_filter_flat(filt_order) takes two 4F flats 
         and gives I-band the logner int quadrant
    For four filter, dome flat was made with two integration times, 
        20(2) for BVR 
        60(45) for I
    """
    
    util.mkdir(calib_dir)
    util.mkdir(dark_dir)
    
    filt_order = "RIVB" # position 1
    flat_num = np.arange(88, 92+1)
    
    print("Filter: ", filt_order)
    
    print("Scanning Flat Frames... ")
    flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(screen_dir, ss) for ss in flat_num]
    scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(screen_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    
    print("Scanning Dark Frames... ")
    dark_num = np.arange(94, 102+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    reduce_STA.treat_overscan(dark_frames)
    
    print("Making Flats... ")
    calib.makeflat(scan_flat_frames, scan_dark_frames[:len(scan_flat_frames)], 
                   f'{calib_dir}flat_{filt_order}.fits', darks=True, fourfilter=True)
    print("Making mask... ")
    calib.make_mask(f'{calib_dir}flat_{filt_order}.fits', f'{calib_dir}mask_{filt_order}.fits',
                       mask_min=0.5, mask_max=1.8,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
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



