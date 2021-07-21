## reduce_2021_05_03.py
## edited by Eden McEwen
## June 2021

import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
from imaka.reduce import reduce_fli
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

root_dir = '/g/lu/data/imaka/onaga/20210503/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Fld2/'
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

# Junk files -- see logs
# Note all in 1x1 binning.

dict_suffix = {'open': '_o',
               'LS': 'LS_c',
               'docz': 'docz2_c',
               'docm': 'docm3_c'}

dict_images = {'open':      [8, 10, 13, 15, 18, 20, 23, 25, 28, 30, 35, 40, 42, 45, 47, 50, 52, 55],
               'LS':        [7, 12, 17, 22, 27, 32, 39, 44, 49, 54],
               'docz':      [9, 14, 19, 24, 29, 34, 41, 46, 51],
               'docm':      [11, 16, 21, 26, 31, 36, 43, 48, 53]}

dict_skies = {'open':      'fld2_sky.fits',
              'LS':        'fld2_sky.fits',
              'docz':      'fld2_sky.fits',
              'docm':      'fld2_sky.fits'}

dict_fwhm = {'open': 14,
             'LS':   6,
             'docz': 6,
             'docm': 6}


def make_flat(): 
    """
    Flats made with twilights
    """
    util.mkdir(calib_dir)
    
    flat_num = np.arange(56, 68+1) # Change per night
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num] # file list
    
    ## STA: needs to deal with overscan
    reduce_STA.treat_overscan(flat_frames) # run once
    
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num] # new file list

    ## Main Function
    calib.makeflat(scan_flat_frames, None, calib_dir + 'flat.fits', darks=False)

    ## Lets also make a mask to use when we call find_stars.
    ## This mask tells us where not to search for stars.
    calib.make_mask(calib_dir + 'flat.fits', calib_dir + 'mask.fits',
                    mask_min=0.7, mask_max=1.6,
                    left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(0, 6+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)

    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky.fits')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
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
        reduce_fli.clean_images(scn_files, out_dir, rebin=1,
                                    sky_frame=sky_dir + sky,
                                    flat_frame=calib_dir + "flat.fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return


def find_stars_fld2():
    ## Loop through all the different data sets
    #for key in dict_suffix.keys():
    #for key in ['docz']:
    for key in ['LS', 'docm']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        reduce_fli.find_stars(img_files, fwhm=fwhm, threshold=6 , N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask.fits')
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)

    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Catalog: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        stats_file = stats_dir + 'stats_' + key + '.fits'
        
        reduce_fli.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)

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
    return None


def stack():
    util.mkdir(stacks_dir)

    ## Loop through all the different data sets
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
        
        reduce_fli.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                               mask_file=calib_dir + 'mask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    reduce_fli.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # reduce_fli.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return