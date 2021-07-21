## reduce_2021_04_30.py
## edited by Eden McEwen
## June 2021

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
massdimm_dir = root_dir + 'reduce/massdimm/'

# skies for bin2 are pulled from previous night
root_dir_2 = '/g/lu/data/imaka/onaga/20210429/sta/'
data_dir_bin2 = root_dir_2 + 'Fld2/'
sky_dir_bin2 = root_dir_2 + 'reduce/sky/'


# This night had two types of binning,
# Suffixes bin1, bin2 has been added to indicate different science image sizes
# DON'T use this night as a template!

dict_suffix = {'open_bin1': '_o',
               'LS_bin1':   'LS_c',
               'docz_bin1': 'docz2_c',
               'open_bin2': '_o',
               'LS_bin2':   'LS_c',
               'docz_bin2': 'docz2_c'}

dict_images = {'open_bin1': [68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94],
               'LS_bin1':   [67, 71, 75, 79, 83, 87, 91],
               'docz_bin1': [69, 73, 77, 81, 85, 89, 93],
               'open_bin2': [64, 66],
               'LS_bin2':   [63],
               'docz_bin2': [65]}

# have bin1 skies only, borrow bin2 from night 0429
dict_skies = {'open_bin1':    'fld2_sky_180_bin1.fits',
              'LS_bin1':      'fld2_sky_180_bin1.fits',
              'docz_bin1':    'fld2_sky_180_bin1.fits',
              'open_bin2':    'fld2_sky_180_bin2.fits',
              'LS_bin2':      'fld2_sky_180_bin2.fits',
              'docz_bin2':    'fld2_sky_180_bin2.fits'}

dict_fwhm = {'open_bin1': 7,
               'LS_bin1':     3,
               'docz_bin1':   3,
               'open_bin2':   7,
               'LS_bin2':     3,
               'docz_bin2':   3}
    

def make_flat(): 
    """
    Makes a flat and mask from flats. 
    Done for two different binnings
    """
    util.mkdir(flat_dir)
    
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

    # bin2 skys were taken 0429, copy files
    shutil.copyfile(sky_dir_bin2+'fld2_sky_180.fits', sky_dir+'fld2_sky_180_bin2.fits')
    shutil.copyfile(sky_dir_bin2+'fld2_sky_180.list', sky_dir+'fld2_sky_180_bin2.list')
    
    return


def reduce_fld2():
    
    util.mkdir(out_dir)

    ## Loop through all datasets
    for key in dict_suffix.keys():   
        bin_num = 'bin1' if 'bin1' in key else 'bin2'  # key should contain bin info
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [data_dir + 'sta{img:03d}{suf:s}.fits'.format(img=ii, suf=suf) for ii in img]
        scn_files = [data_dir + 'sta{img:03d}{suf:s}_scan.fits'.format(img=ii, suf=suf) for ii in img]
        
        reduce_STA.treat_overscan(img_files)
        redu.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=flat_dir+"flat_"+bin_num+".fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return


def find_stars_fld2():
    
    ## Loop through all datasets
    for key in dict_suffix.keys():
    #for key in ['docz']:
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        bin = 'bin1' if 'bin1' in key else 'bin2'  # key should contain bin info 

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        redu.find_stars(img_files, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask_'+bin+'.fits')
        #reduce_fli.find_stars(img_files, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
        #                      mask_flat=flat_dir+"flat"+bin+".fits", mask_min=0.8, mask_max=1.4,
        #                      left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)

                          
    return


def calc_star_stats():
    util.mkdir(stats_dir)
    
    ## Loop through all datasets
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Catalog: ', img)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        stats_file = stats_dir + 'stats_' + key + '.fits'
        redu.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)

    # DEBUG - single threaded
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
            massdimm.append_mass_dimm(stats, massdimm_dir)
        else:
            print('Skipping ' + stats)

    return


def stack():
    util.mkdir(stacks_dir)
    
    ## Loop through all datasets
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
    all_images = []
    
    ## Loop through all datasets
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf + '.fits']
        all_images.append(image_file[0])
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir + 'mask.fits')

    # Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    # DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        
    return
