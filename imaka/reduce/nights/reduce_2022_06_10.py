## reduce_2022_06_10.py
##########################
## edited by Eden McEwen
## June 2022
## An I band run only
## This night had good seeing

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
import re
from imaka.reduce import massdimm
import matplotlib
# matplotlib.use('Agg')

night = '20220610'
root_dir = f'/g/lu/data/imaka/onaga/{night}/sta/'

data_dir = root_dir + 'Fld2/'
twi_dir = root_dir + 'twilights/'
cal_dir = root_dir + 'Calunit/'
dome_dir = root_dir + 'Dome/'
dark_dir = root_dir + 'dark/'

sky_dir = root_dir + 'reduce/sky/' 
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
massdimm_dir = root_dir + 'reduce/massdimm/'

## Junk files -- see logs
## 130 - oversaturated DM voltages

dict_suffix = {'open':  '_o',
               'LS_4':  'LS4WFS_c',
               'LS_3S': 'LS3WFSS_c',
               'LS_3L': 'LS3WFSL_c',
               'TT_3':  'TT3WFSL_c',
              }

dict_images = {'open':  [3,9,15,21,27,33,39,48,51,54,57,69,63,66,69,72],
               'LS_4':  [4,10,16,22,28,34,40],
               'LS_3S': [2,8,14,20,26,32,38],
               'LS_3L': [0,6,12,18,24,30,36,47,50,53,56,59,62,65,68,71],
               'TT_3':  [1,7,13,19,25,31,37,49,52,55,58,61,64,67,70,73],
              }

dict_fwhm = {'open':    15,
               'LS_4':  7,
               'LS_3S': 7,
               'LS_3L': 7,
               'TT_3':  10,
              }


###############################################
### REDUCTION
###############################################

def make_dark():

    util.mkdir(calib_dir)

    ## COPY FROM OLD NIGHT
    shutil.copyfile(root_dir + '../../20220609/sta/reduce/calib/fld2_dark120.fits', calib_dir + 'fld2_dark120.fits')
    
    ## CREATING A DARK
    #print("I dark (120)")
    #dark_num = np.arange(110, 114+1)
    #dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    #scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    
    #reduce_STA.treat_overscan(dark_frames)
    #calib.makedark(scan_dark_frames, calib_dir + 'fld2_dark120.fits')
    
    return


def make_flat(): 
    """
    Makes flat and data mask. 
    Just for single filter I band
    """
    util.mkdir(calib_dir)
    
    ## Copy flat from a previous night
    ## COPY FROM OLD NIGHT
    shutil.copyfile(root_dir + '../../20220609/sta/reduce/calib/domeflat.fits', calib_dir + 'domeflat.fits')
    #shutil.copyfile(root_dir + '../../20210724/sta/reduce/calib/flat_BRIV.fits', calib_dir + 'flat_BRIV.fits')
    
    ## Creating flat from range, I band only
    #flat_num = np.arange(105, 109+1)
    #flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(dome_dir, ss) for ss in flat_num]
    #scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    #reduce_STA.treat_overscan(flat_frames)
    #For darks = use master dark as many times as there are flat frames => reduces noise
    #scan_dark_frames = [calib_dir + 'fld2_dark120.fits' for ss in flat_num]
    #calib.makeflat(scan_flat_frames, scan_dark_frames, calib_dir + 'domeflat.fits', darks=True)
    
    calib.make_mask(calib_dir + 'domeflat.fits', calib_dir + 'domemask.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    return


def make_sky():

    util.mkdir(sky_dir)

    ## COPY A SKY => use a dark
    # shutil.copyfile(root_dir + '../../20210724/sta/dark/dark_044_scan.fits', sky_dir + 'fld2_sky_tmp.fits')
    
    ## CREATING A SKY
    sky_num = np.arange(42,46+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky.fits')

    return


def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    #for key in ['open_RIVB', 'LS_RIVB', 'docz_RIVB']: ## Single key setup
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        redu.clean_images(scn_files, out_dir, rebin=1,
                                    sky_frame=sky_dir + f"fld2_sky.fits",
                                    flat_frame=calib_dir + f"domeflat.fits")#,
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
        fwhm = dict_fwhm[key]
        
        thrsh = 10
        peak_max = 30000
        sharp_lim = 0.8

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]

        redu.find_stars(img_files, fwhm=fwhm, threshold=thrsh, N_passes=2, plot_psf_compare=False,
                        mask_file=calib_dir+f'domemask.fits', peak_max=peak_max, sharp_lim=sharp_lim)
        
    ## DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)
    
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    #for key in ['open_IVBR', 'LS_IVBR', 'docz_IVBR']:
    #for key in []:
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
    #key_i = 'open_BRIV'
    #fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    #image_file = fmt.format(dir=out_dir, img=dict_images[key_i][0], suf=dict_suffix[key_i]) 
    #stats_file = f'{stats_dir}stats_{key_i}.fits'
    #redu.calc_star_stats([image_file], output_stats=stats_file)
    #moffat.fit_moffat_single(image_file,image_file.replace('.fits', '_stars_stats.fits'), 0.2)

        
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
    ## EDITED FOR 4F DATA

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
    ## EDITED FOR 4F DATA
    
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    all_images = []
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]
        
        thrsh = 10
        peak_max = 30000
        sharp_lim = 0.9

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf +'.fits']
        all_images.append(image_file[0])
        
        #redu.find_stars(image_file, fwhm=fwhm, threshold=6, N_passes=2, plot_psf_compare=False, mask_file=calib_dir + f'mask_{filt}.fits')
        redu.find_stars(image_file, fwhm=fwhm, threshold=thrsh, N_passes=2, plot_psf_compare=False,
                        mask_file=calib_dir+f'domemask.fits', peak_max=peak_max, sharp_lim=sharp_lim)

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return
            