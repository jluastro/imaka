 ## reduce_2022_01_25.py
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

night = '20220125'
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


## Junk files - see logs
dict_suffix = {'LS_3wfs_s': 'n3wfs_c',
               'LS_3wfs_w': 'n3wide_c',
               'LS_5wfs': 'n5wfs_c',
               'open':    '_o',
              }

dict_images = {'LS_5wfs':   [0,4,8,10,14,18,22,26,30,34],
               'LS_3wfs_s': [1,5,11,15,19,23,27,31,35],
               'LS_3wfs_w': [2,6,12,16,20,24,28,32,36],
               'open':      [3,7,13,17,21,25,29,33,37],
              }

dict_fwhm = {'LS_5wfs':    8,
               'LS_3wfs_s':  8,
               'LS_3wfs_w':   8,
               'open':      15
            }

###############################################
### REDUCTION
###############################################


def make_flat(): 
    """
    Makes flat and data mask. 
    Just for single filter I band
    CHANGES: took domeflats 01_22
    """
    util.mkdir(calib_dir)
    
    ## Copy flat from previous night because twilight exposures, did not bring flat lamp
    # opted to not do this, no other flat had wfs in the right position
    shutil.copyfile(root_dir + '../../20220122/sta/reduce/calib/domeflat.fits', calib_dir + 'domeflat.fits')
    
    # Mask for masking find_stars areas
    calib.make_mask(calib_dir + 'domeflat.fits', calib_dir + 'domemask.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(38, 44+1) #sky set 1 middle of the night
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir + 'beehive_sky.fits')
    
    return


def reduce_beehive():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
    #for key in ['open']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = 'beehive_sky.fits'

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)    
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        #reduce_STA.treat_overscan_working(img_files)  #BUG
        redu.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=calib_dir + "domeflat.fits")
    return

def rebin_mask():
    # rebin flat
    redu.write_rebin_single(calib_dir+'domeflat.fits', 2, calib_dir+'domeflat_bin2.fits')
    # Mask from rebinned flat
    calib.make_mask(calib_dir + 'domeflat_bin2.fits', calib_dir + 'domemask_bin2.fits',
                       mask_min=0.8*4, mask_max=1.4*4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
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

###############################################
### ANALYSIS
###############################################

def find_stars_beehive():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
    #for key in ['open']:

        img = dict_images[key]
        suf = dict_suffix[key]
        sky = sky_dir + 'beehive_sky.fits'
        
        # o/c loop distinction
        fwhm = 15 if re.search('open', key) else 8
        thrsh = 6 if re.search('open', key) else 7

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)

        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        # Taken from working branch version args
        redu.find_stars(img_files, fwhm=fwhm,  threshold = thrsh, plot_psf_compare=False, mask_file=calib_dir+'domemask.fits')
        
    ## DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)
    
    ## Loop through all the different data sets
    for key in dict_suffix.keys():
    #for key in ['open']:
        
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


def stack_beehive():

    util.mkdir(stacks_dir)

    # Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
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
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        
        # o/c loop distinction
        fwhm = 15 if re.search('open', key) else 8
        thrsh = 6 if re.search('open', key) else 7

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))
        print('    thrsh: ', str(thrsh))

        image_file = [stacks_dir + 'beehive_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=thrsh, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir + 'domemask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return


