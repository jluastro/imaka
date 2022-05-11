 ## reduce_2022_01_22.py
##########################
## edited by Eden McEwen
## January 2022
## Single filter winter cluster
## pulled reduction code from 08_2021

import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
#import scipy
import glob
import re
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

night = '20220122'
root_dir = f'/g/lu/data/imaka/onaga/{night}/sta/'

data_dir = root_dir + 'Beehive-W/'
twi_dir = root_dir + 'twilights/'
cal_dir = root_dir + 'calunit/'
dome_dir = root_dir + 'dome/'
dark_dir = root_dir + 'dark/'
screen_dir = root_dir + 'Screen/'

sky_dir = root_dir + 'reduce/sky/' 
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Beehive-W/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
massdimm_dir = root_dir + 'reduce/massdimm/'

rebin=True

## Junk files - 15 - see logs
dict_suffix = {'LS_3wfs_s_1': 'n3wfs_c',
               'LS_3wfs_w_1': 'n3wide_c',
               'LS_5wfs_1': 'n5wfs_c',
               'open_1':    '_o',
               'LS_3wfs_s_2': 'n3wfs_c',
               'LS_3wfs_w_2': 'n3wide_c',
               'LS_5wfs_2': 'n5wfs_c',
               'open_2':    '_o',
              }

dict_images = {'LS_5wfs_1':   [11,19,23,27,31,35,39,43,47,51,55,59],
               'LS_3wfs_s_1': [12,16,20,24,28,32,36,40,44,48,52,56,60],
               'LS_3wfs_w_1': [13,17,21,25,29,33,37,41,45,49,53,57,61],
               'open_1':      [14,18,22,26,30,34,38,42,46,50,54,58,62],
               'LS_5wfs_2':   [70,74,78,82,86,90,94,98,102,106],
               'LS_3wfs_s_2': [71,75,79,83,87,91,95,99,103,107],
               'LS_3wfs_w_2': [72,76,80,84,88,92,96,100,104,108],
               'open_2':      [73,77,81,85,89,93,97,101,105,109],
              }

dict_sky = {'LS_3wfs_s_1':    'beehive_sky1.fits',
               'LS_3wfs_w_1': 'beehive_sky1.fits',
               'LS_5wfs_1':   'beehive_sky1.fits',
               'open_1':      'beehive_sky1.fits',
               'LS_3wfs_s_2': 'beehive_sky2.fits',
               'LS_3wfs_w_2': 'beehive_sky2.fits',
               'LS_5wfs_2':   'beehive_sky2.fits',
               'open_2':      'beehive_sky2.fits',
              }

dict_fwhm = {'LS_3wfs_s_1':    10,
               'LS_3wfs_w_1':  10,
               'LS_5wfs_1':   10,
               'open_1':      20,
               'LS_3wfs_s_2': 10,
               'LS_3wfs_w_2': 10,
               'LS_5wfs_2':   10,
               'open_2':      20,
              }

###############################################
### REDUCTION
###############################################


def make_flat(): 
    """
    Makes flat and data mask. 
    Just for single filter I band
    CHANGES: took domeflats this night
    """
    util.mkdir(calib_dir)
    
    # Flats: at 200 seconds
    flat_num = np.arange(1, 5+1)
    flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(dome_dir, ss) for ss in flat_num]
    scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    
    # Darks: at 210 seconds
    dark_num = np.arange(6, 10+1)
    dark_frames = ['{0:s}sta{1:03d}_o.fits'.format(dome_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(dome_dir, ss) for ss in dark_num]
    reduce_STA.treat_overscan(dark_frames)
    
    # Make flat
    calib.makeflat(scan_flat_frames, scan_dark_frames, calib_dir + 'domeflat.fits', darks=True)
    
    # Mask for masking find_stars areas
    calib.make_mask(calib_dir + 'domeflat.fits', calib_dir + 'domemask.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)

def make_sky():

    util.mkdir(sky_dir)

    #sky_num = np.arange(63, 69+1) #sky set 1 middle of the night
    #sky_num = np.arange(110, 116+1) #sky set 2 end of the night
    #sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    #reduce_STA.treat_overscan(sky_frames)
    #scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    #calib.makedark(scan_sky_frames, sky_dir + 'beehive_sky2.fits')
    
    # ALl skies
    sky_num = np.append(np.arange(63, 69+1), [110, 111, 113, 114, 115, 116]) # fill skyset
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir + 'beehive_sky_no112.fits')
    
    return


def reduce_beehive():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
    #for key in ['open']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_sky[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)    
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        #reduce_STA.treat_overscan_working(img_files)  #BUG
        redu.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=calib_dir + "domeflat.fits")
    return

def rebin_data():
    """
    Takes the reduced files, and rebins them
    this changes the sape of the file and the header. 
    """
    util.mkdir(out_dir + "bin2/")
    binfac = 2
    for key in dict_suffix.keys():
    #for key in ["LS_3wfs_s_1"]:
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Rebinning: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        rebin_files = [out_dir + 'bin2/sta{img:03d}{suf:s}_scan_clean_bin2.fits'.format(img=ii, suf=suf) for ii in img]
        
        redu.write_rebin(img_files, rebin_files, binfac)
    return

def rebin_mask():
    # rebin flat
    redu.write_rebin_single(calib_dir+'domeflat.fits', 2, calib_dir+'domeflat_bin2.fits')
    # Mask from rebinned flat
    calib.make_mask(calib_dir + 'domeflat_bin2.fits', calib_dir + 'domemask_bin2.fits',
                       mask_min=0.8*4, mask_max=1.4*4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    return
    

###############################################
### ANALYSIS
###############################################

def find_stars_beehive():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
    #for key in ['LS_3wfs_s_1']:

        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_sky[key]
        
        fwhm = dict_fwhm[key]
        sharp_lim = 0.9 #dict_sharp[key]
        thrsh = 4
        peak_max = 30000

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        print('    Rebin: ', rebin_data)
        
        if rebin_data:
            img_files = [out_dir + 'bin2/sta{img:03d}{suf:s}_scan_clean_bin2.fits'.format(img=ii, suf=suf) for ii in img]
            mask_f = calib_dir+'domemask_bin2.fits'
            fwhm = fwhm/2
        else:
            img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
            mask_f = calib_dir+'domemask.fits'
            
        # find stars on a starlist in parallel
        redu.find_stars(img_files, fwhm=fwhm,  threshold = thrsh, plot_psf_compare=False, mask_file=mask_f, peak_max=peak_max, sharp_lim=sharp_lim)
        
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
        
        if rebin_data:
            img_files = [out_dir + 'bin2/sta{img:03d}{suf:s}_scan_clean_bin2.fits'.format(img=ii, suf=suf) for ii in img]
            stats_file = stats_dir + 'stats_' + key + '_bin2.fits'
        else:
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


def stack_beehive():

    util.mkdir(stacks_dir)

    # Loop through all the different data sets and reduce them.
    #for key in dict_suffix.keys():
    for key in ['open_2', 'LS_3wfs_s_2', 'LS_3wfs_w_2', 'LS_5wfs_2']:
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        output_root = stacks_dir + 'beehive_stack_' + suf
        redu.shift_and_add(img_files, starlists, output_root, method='mean')
        
    return


def analyze_stacks():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    all_images = []
    #for key in dict_suffix.keys():
    for key in ['open_2', 'LS_3wfs_s_2', 'LS_3wfs_w_2', 'LS_5wfs_2']:
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'beehive_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir + 'domemask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return

