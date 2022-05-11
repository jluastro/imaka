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
import re
from imaka.reduce import massdimm
from imaka.reduce import reduce_STA
import matplotlib
matplotlib.use('Agg')

night = '20200121'
targt = 'Beehive-W'

root_dir = f'/g/lu/data/imaka/onaga/{night}/sta/'
sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + f'{targt}/'
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + f'reduce/{targt}/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

# pulling calib files from following night
root_alt_dir = '/g/lu/data/imaka/onaga/20200122/sta/'
data_alt_dir = root_alt_dir + f'{targt}/'
twi_alt_dir = root_alt_dir + 'twilights/'

fnum_o_I30 = [38, 40, 42, 44, 46]
fnum_o_I60 = [48, 50, 52, 54, 56]
fnum_o_4F_1 = [65, 67, 70, 79, 81]
fnum_o_4F_2 = [83, 85, 87, 89, 91, 93]
fnum_o_4F_3 = [107, 109, 111, 113, 115]

fnum_c_I30 = [36, 37, 39, 41, 43, 45]
fnum_c_I60 = [47, 49, 51, 53, 55]
fnum_c_4F_1 = [64, 66, 69, 78, 80]
fnum_c_4F_2 = [82, 84, 86, 88, 90, 92]
fnum_c_4F_3 = [106, 108, 110, 112, 114]

dict_suffix = {'o_60': '_o',
               'o_30': '_o',
               'o_4F_1': '_o',
               'o_4F_2': '_o',
               'o_4F_3': '_o',
               'c_60': 'x10LS5WFS_c',
               'c_30': 'x10LS5WFS_c',
               'c_4F_1': 'x10LS5WFS_c',
               'c_4F_2': 'x10LS5WFS_c',
               'c_4F_3': 'x10LS5WFS_c'}
    
dict_images = {'o_60': fnum_o_I30,
               'o_30': fnum_o_I60,
               'o_4F_1': fnum_o_4F_1,
               'o_4F_2': fnum_o_4F_2,
               'o_4F_3': fnum_o_4F_3,
               'c_60': fnum_c_I30,
               'c_30': fnum_c_I60,
               'c_4F_1': fnum_c_4F_1,
               'c_4F_2': fnum_c_4F_2,
               'c_4F_3': fnum_c_4F_3}

dict_skies =  {'o_60': 'beehive_sky_60.fits',
               'o_30': 'beehive_sky_60.fits', # BUG
               'o_4F_1': 'beehive_sky_60.fits',
               'o_4F_2': 'beehive_sky_60.fits',
               'o_4F_3': 'beehive_sky_60.fits',
               'c_60': 'beehive_sky_60.fits',
               'c_30': 'beehive_sky_60.fits', # BUG
               'c_4F_1': 'beehive_sky_60.fits',
               'c_4F_2': 'beehive_sky_60.fits',
               'c_4F_3': 'beehive_sky_60.fits'}


def make_flat(): 

    util.mkdir(flat_dir)

    # We don't have one for this run, so I just symbolically linked from an really
    # really old run.
    # In sta/reduce/calib/
    # ln -s /g/lu/data/imaka/onaga/20181223/sta/reduce/calib/flat.fits
    
    # 2021 update: will link to flats in 20200122 (the following night)
    flat_num = np.arange(221, 230+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_alt_dir, ss) for ss in flat_num]
    reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_alt_dir, ss) for ss in flat_num]
    calib.makeflat(scan_flat_frames, [], calib_dir + 'flat.fits', darks=False)
    
    # Mask for masking find_stars areas
    calib.make_mask(calib_dir + 'flat.fits', calib_dir + 'mask.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)

    return


def make_sky():
    
    util.mkdir(sky_dir)
 
    # 60s exposure
    sky_num = np.arange(208, 217+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_alt_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_alt_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'beehive_sky_60.fits')
    
    # 150s exposure
    sky_num = np.arange(218, 220+1) 
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_alt_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)
    scan_sky_frames = ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_alt_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir+'beehive_sky_150.fits')
    
    return


def reduce_beehive():

    util.mkdir(out_dir)

    ## Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
    #for key in ['c_30']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)    
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        #reduce_STA.treat_overscan(img_files)
        reduce_STA.treat_overscan_working(img_files)  #BUG
        reduce_fli.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=calib_dir + "flat.fits")
    return


def find_stars_beehive():
    
    # Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():

        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        
        # o/c loop distinction
        fwhm = 5 if re.search('c_', key) else 8
        thrsh = 10 if re.search('c_', key) else 5

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)

        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        # Taken from working branch version args
        reduce_fli.find_stars(img_files, fwhm=fwhm,  threshold = thrsh, plot_psf_compare=False, mask_file=calib_dir+'mask.fits')

    # DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
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

    # DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0])
    # stats_file = stats_dir + 'stats_LS_c.fits'
    # redu.calc_star_stats(image_file, stats_file, flux_percent=0.2)

    return


def append_massdimm():

    massdimm.fetch_data('20200121', massdimm_dir) ## not working, http doesn't exist
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

    # Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        starlists = [out_dir + 'sta{img:03d}{suf:s}_scan_clean_stars.txt'.format(img=ii, suf=suf) for ii in img]
        output_root = stacks_dir + 'beehive_stack_' + suf
        reduce_fli.shift_and_add(img_files, starlists, output_root, method='mean')
        
    return


def analyze_stacks():
    # Loop through all the different data sets and reduce them.
    all_images = []
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'beehive_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        reduce_fli.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir + 'mask.fits')

    # Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    reduce_fli.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    # DEBUG - single threaded
    # image_file = stacks_dir + 'beehive_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

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

## Splitting the initial dictonary keys above means they don't need to be again split here.
## Intentionally redundant here

filters = ['B', 'V', 'R', 'I']

# Open
rot_1_o = [65, 67, 70, 79, 81]      # key: o_4F_1
rot_2_o = [83, 85, 87, 89, 91, 93]  # key: o_4F_2
rot_3_o = [107, 109, 111, 113, 115] # key: o_4F_3

# Closed
rot_1_c = [61, 64, 66, 69, 78, 80]  # key: c_4F_1
rot_2_c = [82, 84, 86, 88, 90, 92]  # key: c_4F_2
rot_3_c = [106, 108, 110, 112, 114] # key: c_4F_3

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
            reduce_fli.calc_star_stats(img_files, output_stats=stats_file, starlists=starlists, fourfilt=True)
            print("Starting moffat fitting")
            moffat.fit_moffat(img_files, stats_file, starlists=starlists)
            
            ## SINGLE THREAD
            # reduce_fli.calc_star_stats_single(img_files[0], starlists[0], True)
    
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
