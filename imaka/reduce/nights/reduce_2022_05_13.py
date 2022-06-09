## reduce_2022_05_13.py
##########################
## edited by Eden McEwen
## May 2022
## A four filter run
## Did position 3 (VBRI) this night, not planning on rotating

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

night = '20220513'
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

dict_suffix = {'open': '_o',
               'LS':   'LS_c',
              }

dict_images = {'LS':   [86,88,90,92,94,96,103,115,107,109,111,113,115,122,124,126,128,130,132,134,136,138,142], 
               'open': [87,89,91,93,95,97,104,106,108,110,112,114,116,123,125,127,129,131,133,135,137,139,143], 
              }

dict_fwhm = {'open': 12,
             'LS': 8,
            }  

# only include filter if key was a 4F file
dict_filt = {'open': 'VBRI',
             'LS':   'VBRI',
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
    #shutil.copyfile(root_dir + '../../20210724/sta/reduce/calib/flat_BRIV.fits', calib_dir + 'flat_BRIV.fits')
    #shutil.copyfile(root_dir + '../../20210723/sta/reduce/calib/flat_RIVB.fits', calib_dir + 'flat_RIVB.fits')
    
    ## Creating flat from range, I band only
    #flat_num = np.arange(37, 49+1)
    #flat_frames = ['{0:s}dome_{1:03d}.fits'.format(dome_dir, ss) for ss in flat_num]
    #scan_flat_frames = ['{0:s}dome_{1:03d}_scan.fits'.format(dome_dir, ss) for ss in flat_num]
    #reduce_STA.treat_overscan(flat_frames)
    #calib.makeflat(scan_flat_frames, None, calib_dir + 'domeflat_I.fits', darks=False)
    
    calib.make_mask(calib_dir + 'flat_VBRI.fits', calib_dir + 'mask_VBRI.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    return


def make_sky():

    util.mkdir(sky_dir)

    ## COPY A SKY => use a dark
    # shutil.copyfile(root_dir + '../../20210724/sta/dark/dark_044_scan.fits', sky_dir + 'fld2_sky_tmp.fits')
    
    ## CREATING A SKY
    sky_num = np.arange(117, 121+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky_VBRI.fits')

    return

def make_dark():

    util.mkdir(sky_dir)

    ## COPY A SKY => use a dark
    # shutil.copyfile(root_dir + '../../20210724/sta/dark/dark_044_scan.fits', sky_dir + 'fld2_sky_tmp.fits')
    
    ## CREATING A SKY
    print("I quad flat (60)")
    dark_num = np.arange(81, 85+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    
    reduce_STA.treat_overscan(dark_frames)
    calib.makedark(scan_dark_frames, calib_dir + 'fld2_dark60_VBRI.fits')
    
    print("BVR quad flat (20)")
    dark_num = np.arange(76,80+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    
    reduce_STA.treat_overscan(dark_frames)
    calib.makedark(scan_dark_frames, calib_dir + 'fld2_dark20_VBRI.fits')
    
    return

def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    #for key in ['open_RIVB', 'LS_RIVB', 'docz_RIVB']: ## Single key setup
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
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]
        filt = dict_filt[key]
        fwhm = dict_fwhm[key]
        
        thrsh = 10
        peak_max = 30000
        sharp_lim = 0.9

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      filt: ', filt)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]

        redu.find_stars(img_files, fwhm=fwhm, threshold=thrsh, N_passes=2, plot_psf_compare=False,
                        mask_file=calib_dir+f'mask_{filt}.fits', peak_max=peak_max, sharp_lim=sharp_lim)
        
    ## DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
    return

def filter_stars():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]
        filt = dict_filt[key]
        fwhm = dict_fwhm[key]

        print('Filtering: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      filt: ', filt)
        
        star_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        # TODO: move this into a faster system, in the reduce file
        for sf in star_files:
            stars = table.Table.read(sf, format='ascii')
            
            formats = {'xcentroid': '%8.3f', 'ycentroid': '%8.3f', 'sharpness': '%.2f',
                   'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
                   'flux': '%10.6f', 'mag': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
                   'theta': '%6.3f'}
            
            stars.write(sf.split(".txt")[0]+"_orig.txt", format='ascii.fixed_width', 
                        delimiter=None, bookend=False, formats=formats, overwrite=True)
            
            stars = stars[stars['peak'] > 100]
            stars = stars[stars['peak'] < 30000]
            stars = stars[stars['sharpness'] < 0.7]

            stars.write(sf, format='ascii.fixed_width', delimiter=None, bookend=False, formats=formats, overwrite=True)
                          
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
        filt = dict_filt[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        output_root = stacks_dir + 'fld2_stack_' + suf + '_' + filt ## EDITED LINE
        
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
        filt = dict_filt[key]
        
        thrsh = 10
        peak_max = 30000
        sharp_lim = 0.9

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf +  '_' + filt +'.fits'] ## EDITED LINE
        all_images.append(image_file[0])
        
        #redu.find_stars(image_file, fwhm=fwhm, threshold=6, N_passes=2, plot_psf_compare=False, mask_file=calib_dir + f'mask_{filt}.fits')
        redu.find_stars(image_file, fwhm=fwhm, threshold=thrsh, N_passes=2, plot_psf_compare=False,
                        mask_file=calib_dir+f'mask_{filt}.fits', peak_max=peak_max, sharp_lim=sharp_lim)

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

filters = ['V', 'R', 'I', 'B']

## Catching any missed name changes
dict_orders_rot = dict_filt
dict_images_rot = dict_images
dict_suffix_rot = dict_suffix


def make_flat_4F(): 
    """
    Makes flat and data mask. 
    For four filter, dome flat was made with two integration times, 20(2) for BVR 60(45) for I
    """
    util.mkdir(calib_dir)
    
    ## Darks are the same all night
    filt_order = "VBRI"
    flat_num_60 = np.arange(63, 67+1)
    flat_num_20 = np.array([73, 74, 75, 73, 74])
    
    print(filt_order)
    
    print("I quad flat (60)")
    dark_num = np.arange(81, 85+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(dome_dir, ss) for ss in flat_num_60]
    scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(dome_dir, ss) for ss in flat_num_60]
    
    reduce_STA.treat_overscan(dark_frames)
    reduce_STA.treat_overscan(flat_frames)
    calib.makeflat(scan_flat_frames, scan_dark_frames, 
                   f'{calib_dir}domeflat_60_{filt_order}.fits', darks=True, fourfilter=True)
    
    print("BVR quad flat (20)")
    dark_num = np.arange(76,80+1)
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(dark_dir, ss) for ss in dark_num]
    scan_dark_frames = ['{0:s}dark_{1:03d}_scan.fits'.format(dark_dir, ss) for ss in dark_num]
    flat_frames = ['{0:s}sta{1:03d}_o.fits'.format(dome_dir, ss) for ss in flat_num_20]
    scan_flat_frames = ['{0:s}sta{1:03d}_o_scan.fits'.format(dome_dir, ss) for ss in flat_num_20]
    
    reduce_STA.treat_overscan(dark_frames)
    reduce_STA.treat_overscan(flat_frames)
    calib.makeflat(scan_flat_frames, scan_dark_frames, 
                   f'{calib_dir}domeflat_20_{filt_order}.fits', darks=True, fourfilter=True)
    
    # Combining two flats based on filter orientation
    print("Combining: I quad flat (60) & BVR quad flat (20)")
    calib.combine_filter_flat(f'{calib_dir}domeflat_60_{filt_order}.fits',
                              f'{calib_dir}domeflat_20_{filt_order}.fits', 
                              f'{calib_dir}flat_{filt_order}.fits', filt_order, flip_180=True)
    
    calib.make_mask(f'{calib_dir}flat_{filt_order}.fits', f'{calib_dir}mask_{filt_order}.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    return


def split_4F_starlists():
    ## Split all starlists by filter, given rotation order
    ## only 4F files should have a dict_filt key entry
    util.mkdir(out_dir+'4F/')
    for key in dict_filt.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        odr = dict_filt[key]
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        reduce_STA.four_filt_split(starlists, odr, flip_180=True)
    return
    

def calc_4F_stats():
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
                starlists += [out_dir + '4F/sta{img:03d}{suf:s}_scan_clean_{f:s}_{odr:s}_stars.txt'.format(img=ii, suf=suf, f=f, odr=odr) for ii in img]
            starlist_stats = [strlst.replace('_stars.txt', '_stars_stats.fits') for strlst in starlists]
            print(f"Calc Star_Stats: {suf} \n Filter: {f}")
            redu.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
            print("Starting moffat fitting")
            moffat.fit_moffat(img_files, stats_file, starlists=starlist_stats)
            
            ## DEBBUG: SINGLE THREAD
            #print("stars: ", starlists[0])
            #print("stats: ", stats_file)
            #redu.calc_star_stats_single(img_files[0], starlists[0], True)
            #moffat.fit_moffat_single(img_files[0], starlist_stats[0], 0.2)
            #break
        #break 
    return

def update_4F_stats():
    # this function updates the "FILTER" column on the 4F stats summaries
    suffixes = list(set(dict_suffix_rot.values()))
    for suf in suffixes:
        # Iterating through filters
        for f in filters:
            s_f = stats_dir + f'stats_{suf}_{f}.fits'
            print("Updating: ", s_f)
            stats = table.Table.read(s_f)
            # from the number of the file we can add in a column for order. this could have been easier at a different step?
            # list of all the image numbers used in this file
            images = [re.search(f"Fld2/sta(.*){suf}", filt).group(1) for filt in stats["Image"]]
            # want to find the key corresponding to these
            f_ord = [next(key for key, value in dict_images.items() if (int(i) in value)).split("_")[1] for i in images]
            # then get the filter order depending on that key
            stats.replace_column('FILTER', table.Column(np.repeat(f, len(stats))))
            if not 'F_ORD'in stats.colnames:
                stats.add_column(table.Column(f_ord, name='F_ORD'))
            if not 'wavelength' in stats.colnames:
                stats.add_column(table.Column(np.repeat(util.get_wavelength(f), len(stats)), name='wavelength'))
            if not 'quad' in stats.colnames:
                stats.add_column(table.Column([util.get_quad(f, odr) for odr in f_ord], name='quad'))
            stats.write(s_f, overwrite=True)
                             
    return

# must stack and find stars for stacks previously.
def split_4F_stacks():
    util.mkdir(stacks_dir+'4F/')
    all_images = []
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        odr = dict_filt[key]

        starlists = [stacks_dir + 'fld2_stack_' + suf +  '_' + odr +'_stars.txt'] # unfortunate naming error for stacks
        reduce_STA.four_filt_split(starlists, odr, flip_180=True)
    return
    

def analyze_4F_stacks():
    ## Loop through all the different data sets
    #for key in ['set_name']: ## Single key setup
    all_images = []
    all_starlists = []
    
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        odr = dict_filt[key]
        
        all_images += [stacks_dir + f'fld2_stack_{suf}_{odr}.fits' for f in filters]
        all_starlists += [stacks_dir + f'4F/fld2_stack_{suf}_{odr}_{f}_{odr}_stars.txt' for f in filters] # unfortunate naming error
    
    stats_file = stats_dir + 'stats_stacks.fits'
    starlist_stats = [strlst.replace('_stars.txt', '_stars_stats.fits') for strlst in all_starlists]
    ## Calc stats on all the stacked images
    redu.calc_star_stats(all_images, output_stats=stats_file, starlists=all_starlists, fourfilt=True)
    print("Starting moffat fitting")
    moffat.fit_moffat(all_images, stats_file, starlists=starlist_stats, flux_percent=0.2)

    return