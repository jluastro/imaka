## reduce_2021_05_02.py
## edited by Eden McEwen
## June 2021

import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
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

root_dir = '/g/lu/data/imaka/onaga/20210502/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Fld2/'
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
massdimm_dir = root_dir + 'reduce/massdimm/'

# Junk files -- see logs
#  -- old docz
# Note all in 1x1 binning.

dict_suffix = {'open': '_o',
               'LS': 'LS_c',
               'docz': 'docz2_c',
               'doczskycl': 'doczskycl_c'}

dict_images = {'open':      [10, 12, 15, 17, 20, 22, 26, 28, 31, 33, 36, 38],
               'LS':        [9,  14, 19, 25, 30, 35],
               'docz':      [11, 16, 21, 27, 32, 37],
               'doczskycl': [13, 18, 23, 29, 34, 39]}

# These don't exist yet... make some temp files with zeros.
dict_skies = {'open':      'fld2_sky.fits',
              'LS':        'fld2_sky.fits',
              'docz':      'fld2_sky.fits',
              'doczskycl': 'fld2_sky.fits'}

dict_fwhm = {'open': 14,
             'LS':   6,
             'docz': 6,
             'doczskycl': 6}

def make_flat(): 
    """
    Make a flat... this will be with dome flats for now. 
    These are junk... replace with twilight flats. 
    """
    util.mkdir(calib_dir)
    
    ## Copy flight from previous night if twilight exposures
    ## are a little saturated.
    #shutil.copyfile(root_dir + '../../20210430/sta/reduce/calib/flat_bin1.fits', calib_dir + 'flat.fits')
    
    ## assuming these aren't overexposed
    flat_num = np.arange(59, 70+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    #reduce_STA.treat_overscan(flat_frames)
    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
  
    calib.makeflat(scan_flat_frames, None, calib_dir + 'flat.fits', darks=False)

    ## Lets also make a mask to use when we call find_stars.
    ## This mask tells us where not to search for stars.
    ## UPDATE: mask_min and mask_max were hand calculated 6/14/2021
    calib.make_mask(calib_dir + 'flat.fits', calib_dir + 'mask.fits',
                       mask_min=0.5, mask_max=1.8,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(1, 8+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)

    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky.fits')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    ## Loop through all the different data sets
    for key in dict_suffix.keys():
    #for key in ['doczskycl']:
        
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


def find_stars_fld2():
    ## Loop through all the different data sets
    for key in dict_suffix.keys():

        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)

        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        redu.find_stars(img_files, fwhm=fwhm, threshold=8, N_passes=2, plot_psf_compare=False,
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
        
        redu.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)

        # stats1 = table.Table.read(stats_file)
        # redu.add_frame_number_column(stats1)
        # stats1.write(stats_file, overwrite=True)

        # stats_file_mdp = stats_file.replace('.fits', '_mdp.fits')
        # stats2 = table.Table.read(stats_file_mdp)
        # redu.add_frame_number_column(stats1)
        # stats2.write(stats_file_mdp, overwrite=True)

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
        
        redu.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir + 'mask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return