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

root_dir = '/g/lu/data/imaka/onaga/20210501/sta/'

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
               'doczskycl': 'doczskycl_c',
               'z10glao': 'z10glao_c',
               'z10scao': 'Z10scao_c'}

dict_images = {'open':      [73, 75, 78, 80, 83, 85, 88, 90, 93, 95, 103, 105, 109],
               'LS':        [72, 77, 82, 87],
               'docz':      [74, 79, 84, 89],
               'doczskycl': [76, 81, 86, 91],
               'z10glao':   [94, 104],
               'z10scao':   [96, 106, 107]
              }

# These don't exist yet... make some temp files with zeros.
dict_skies = {'open':      'fld2_sky.fits',
              'LS':        'fld2_sky.fits',
              'docz':      'fld2_sky.fits',
<<<<<<< HEAD
              'doczskycl': 'fld2_sky.fits'}

# all bin 1
dict_fwhm = {'open': 14,
             'LS':   6,
             'docz': 6,
             'doczskycl': 6}
=======
              'doczskycl': 'fld2_sky.fits',
              'z10glao': 'fld2_sky.fits',
              'z10scao': 'fld2_sky.fits'
              }

dict_fwhm = {'open': 12,
             'LS': 5,
             'docz': 5,
             'doczskycl': 5,
             'z10glao': 5,
             'z10scao': 5
            }
>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
    

def make_flat(): 
    """
    Make a flat... this will be with dome flats for now. 
    These are junk... replace with twilight flats. 
    """
    util.mkdir(calib_dir)
    
    # Copy flight from previous night because twilight exposures
    # were a little saturated.
    shutil.copyfile(root_dir + '../../20210430/sta/reduce/calib/flat_bin1.fits', calib_dir + 'flat.fits')

    # Lets also make a mask to use when we call find_stars.
    # This mask tells us where not to search for stars.
    # UPDATE: mask_min and mask_max were hand calculated 6/14/2021
    calib.make_mask(calib_dir + 'flat.fits', calib_dir + 'mask.fits',
                       mask_min=0.5, mask_max=1.8,
                       left_slice=20, right_slice=20, top_slice=25, bottom_slice=25)
    
    return


def make_sky():

    util.mkdir(sky_dir)

    sky_num = np.arange(97, 101+1)
    sky_frames = ['{0:s}sky_{1:03d}_o.fits'.format(data_dir, ss) for ss in sky_num]
    reduce_STA.treat_overscan(sky_frames)

    scan_sky_frames =  ['{0:s}sky_{1:03d}_o_scan.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(scan_sky_frames, sky_dir + 'fld2_sky.fits')
    
    return


def reduce_fld2():

    util.mkdir(out_dir)

    # Loop through all the different data sets and reduce them.
    # for key in ['doczskycl']:
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
        redu.clean_images(scn_files, out_dir, rebin=1,
                                    sky_frame=sky_dir + sky,
                                    flat_frame=calib_dir + "flat.fits")#,
                                # fix_bad_pixels=True, worry_about_edges=True)

    return


def find_stars_fld2():
<<<<<<< HEAD
    ## Loop through all the different data sets and reduce them.
=======
    # Loop through all the different data sets and reduce them.
>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
    # for key in ['docz']:
    for key in dict_suffix.keys():
        
        img = dict_images[key]
        suf = dict_suffix[key]
        sky = dict_skies[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('      Sky: ', sky)
        
        img_files = [out_dir + 'sta{img:03d}{suf:s}_scan_clean.fits'.format(img=ii, suf=suf) for ii in img]
<<<<<<< HEAD
        reduce_fli.find_stars(img_files, fwhm=fwhm, threshold=6, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask.fits')
        
=======
        redu.find_stars(img_files, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
                              mask_file=calib_dir+'mask.fits')
        
    # DEBUG - single threaded
    # fmt = '{dir}sta{img:03d}{suf:s}_scan_clean.fits'
    # image_file = fmt.format(dir=out_dir, img=dict_images['LS_c'][0], suf=dict_suffix['LS_c'][0]) 
    # redu.find_stars_single(image_file, dict_fwhm['LS_c'], 3, 2, False, calib_dir + 'mask.fits')
                          
>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
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
<<<<<<< HEAD
        reduce_fli.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)

        # stats1 = table.Table.read(stats_file)
        # reduce_fli.add_frame_number_column(stats1)
        # stats1.write(stats_file, overwrite=True)

        # stats_file_mdp = stats_file.replace('.fits', '_mdp.fits')
        # stats2 = table.Table.read(stats_file_mdp)
        # reduce_fli.add_frame_number_column(stats1)
        # stats2.write(stats_file_mdp, overwrite=True)
=======
        reduce_STA.calc_star_stats(img_files, output_stats=stats_file)
        moffat.fit_moffat(img_files, stats_file, flux_percent=0.2)
>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663

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
<<<<<<< HEAD
    
=======

>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
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
<<<<<<< HEAD
    ## Loop through all the different data sets and reduce them.
=======
    # Loop through all the different data sets and reduce them.
>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
    all_images = []
    for key in dict_suffix.keys():
        img = dict_images[key]
        suf = dict_suffix[key]
        fwhm = dict_fwhm[key]

        print('Working on: {1:s}  {0:s}'.format(key, suf))
        print('   Images: ', img)
        print('     Fwhm: ', str(fwhm))

        image_file = [stacks_dir + 'fld2_stack_' + suf + '.fits']
        all_images.append(image_file)
        
<<<<<<< HEAD
        # reduce_fli.find_stars(image_file, fwhm=fwhm, threshold=3, N_passes=2, plot_psf_compare=False,
        #                       mask_file=calib_dir + 'mask.fits')

    ## Calc stats on all the stacked images
    out_stats_file = stats_dir + 'stats_stacks.fits'
    # reduce_fli.calc_star_stats(all_images, output_stats=out_stats_file)
    moffat.fit_moffat(all_images, out_stats_file, flux_percent=0.2)

    ## DEBUG - single threaded
    # image_file = stacks_dir + 'fld2_stack_' + dict_suffix['open'] + '.fits'
    # reduce_fli.find_stars_single(image_file, dict_fwhm['open'], 3, 2, False, calib_dir + 'mask.fits')        

    return
=======
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

>>>>>>> ccfaf181eb56310de080e8e0a2249198e6823663
