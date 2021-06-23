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

root_dir = '/g/lu/data/imaka/onaga/20210429/sta/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'Fld2/'
calib_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/Fld2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
twi_dir = root_dir + 'twilights/'
old_cals_dir = root_dir + 'cals_from_20210403/'
massdimm_dir = root_dir + 'reduce/massdimm/'

# Junk files -- see logs
# 62, 60 -- old docz
# cmat 0: LS same all night
# cmat 1: zonal (changed after frame 59)
# cmat 2: modal

dict_suffix = {'open_30': '_o',
               'LS_30':   'LS_c',                   
               'docz_30': 'docz2_c',                   
               'open': '_o',
               'LS': 'LS_c',
               'docz': 'docz2_c'}

dict_images = {'open_30': [58, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105],
               'LS_30':   [56, 86, 90, 94, 98, 102],
               'docz_30': [88, 92, 96, 100, 104],
               'open':    [61, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85],
               'LS':      [59, 63, 66, 70, 74, 78, 82],
               'docz':    [64, 68, 72, 76, 80, 84]}

# These don't exist yet... make some temp files with zeros.
dict_skies = {'open_30': 'fld2_sky_30.fits',
              'LS_30':   'fld2_sky_30.fits',
              'docz_30': 'fld2_sky_30.fits',
              'open':    'fld2_sky_180.fits',
              'LS':      'fld2_sky_180.fits',
              'docz':    'fld2_sky_180.fits'}

dict_fwhm = {'open_30': 7,
             'LS_30':   3,
             'docz_30': 3,
             'open': 7,
             'LS': 3,
             'docz': 3}

def make_flat(): 
    """
    Make a flat... this will be with dome flats for now. 
    These are junk... replace with twilight flats. 
    """
    util.mkdir(calib_dir)

    flat_num = np.arange(115, 120+1)
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(twi_dir, ss) for ss in flat_num]
    # reduce_STA.treat_overscan(flat_frames)

    scan_flat_frames = ['{0:s}twi_{1:03d}_scan.fits'.format(twi_dir, ss) for ss in flat_num]
    # calib.makeflat(scan_flat_frames, None, calib_dir + 'flat.fits', darks=False)

    # Lets also make a mask to use when we call find_stars.
    # This mask tells us where not to search for stars.
    calib.make_mask(calib_dir + 'flat.fits', calib_dir + 'mask.fits',
                       mask_min=0.8, mask_max=1.4,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
        
    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num_180 = np.arange(106, 108+1)
    sky_frames_180 = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num_180]
    reduce_STA.treat_overscan(sky_frames_180)

    sky_num_30 = np.arange(109, 111+1)
    sky_frames_30 = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num_30]
    reduce_STA.treat_overscan(sky_frames_30)

    # Put all the 30 sec and 180 sec together. We have no choice since we had so few. 
    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num_30]
    scan_sky_frames += ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num_180]

    rescale_to_180 = [6 for ss in sky_num_30]
    rescale_to_180 += [1 for ss in sky_num_180]
    rescale_to_30 = [1 for ss in sky_num_30]
    rescale_to_30 += [(1./6.) for ss in sky_num_180]

    # We also need to account for the increasing sky brightness.
    for ff in range(len(scan_sky_frames)):
        img = fits.getdata(scan_sky_frames[ff])

        # Fix 180 scales
        img_180 = img * rescale_to_180[ff]
        mean_180, median_180, stddev_180 = sigma_clipped_stats(img_180, sigma_lower=4, sigma_upper=2, maxiters=10)
        if ff == 0:
            mean0_180 = mean_180
        rescale_to_180[ff] *= mean0_180 / mean_180

        # Fix 30 scales
        img_30 = img * rescale_to_30[ff]
        mean_30, median_30, stddev_30 = sigma_clipped_stats(img_30, sigma_lower=4, sigma_upper=2, maxiters=10)
        if ff == 0:
            mean0_30 = mean_30
        rescale_to_30[ff] *= mean0_30 / mean_30
        

    calib.makedark(scan_sky_frames, sky_dir+'fld2_sky_180.fits', rescale=rescale_to_180)
    calib.makedark(scan_sky_frames, sky_dir+'fld2_sky_30.fits', rescale=rescale_to_30)
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    # Loop through all the different data sets and reduce them.
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
        redu.clean_images(scn_files, out_dir, rebin=1, sky_frame=sky_dir + sky, flat_frame=calib_dir + "flat.fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return


def find_stars_fld2():
    # Loop through all the different data sets and reduce them.
    for key in dict_suffix.keys():

        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)

        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
        redu.find_stars(img_files, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask.fits')

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

    # Loop through all the different data sets and reduce them.
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
    # Loop through all the different data sets and reduce them.
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

    # Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    redu.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    # DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # redu.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return

