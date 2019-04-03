import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob, shutil
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
import pdb
import os

imaka_dir = '/Users/jlu/data/imaka/'
root_dir = imaka_dir + '20170218/fli/'

def make_flat():
    # Didn't take flats, so just copy over the one from Friday night.
    # Note this is an R-band flat; but we will use it for R and I.
    old_flat = imaka_dir + '20170214/fli/reduce/calib/flat.fits'
    new_flat = root_dir + 'reduce/calib/flat.fits'

    util.mkdir(root_dir + 'reduce/calib/')
    
    shutil.copyfile(old_flat, new_flat)
    return

def make_sky():
    data_dir = root_dir + 'Pleiades/'
    sky_dir = root_dir + 'reduce/sky/'

    util.mkdir(sky_dir)
    
    # 15 sec integration time
    sky_num = np.arange(73, 83)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky_r_15s.fits')

    # 15 sec integration time
    sky_num = np.arange(130, 146)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky_i_15s.fits')
    
    return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    util.mkdir(out_dir)

    ##########
    # R-band
    ##########
    
    # Open Loop 15 s
    fnum = [31, 38, 45, 52, 59, 66]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame=flat_dir + 'flat.fits')
    
    # Closed Loop - a, 15s
    fnum = [30, 37, 44, 51, 58, 65, 72]
    img_files = [data_dir + 'obj{0:04d}_ca.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - b, 15s
    fnum = [28, 35, 42, 49, 56, 63, 70]
    img_files = [data_dir + "obj{0:04d}_cb.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - c, 15s
    fnum = [27, 34, 41, 48, 55, 62, 69]
    img_files = [data_dir + "obj{0:04d}_cc.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - d, 15s
    fnum = [26, 33, 40, 47, 54, 61, 68]
    img_files = [data_dir + "obj{0:04d}_cd.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - e, 15s
    fnum = [25, 32, 39, 46, 53, 60, 67]
    img_files = [data_dir + "obj{0:04d}_ce.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')
    
    # TT, 15s
    fnum = [29, 36, 43, 50, 57, 64, 71]
    img_files = [data_dir + "obj{0:04d}_tt.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')


    ##########
    # I-band
    ##########
    # Open Loop 15 s
    fnum = [94, 101, 108, 115, 129]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame=flat_dir + 'flat.fits')
    
    # Closed Loop - a, 15s
    fnum = [93, 100, 107, 114, 128]
    img_files = [data_dir + 'obj{0:04d}_ca.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - b, 15s
    fnum = [91, 98, 105, 112, 119]
    img_files = [data_dir + "obj{0:04d}_cb.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - c, 15s
    fnum = [90, 97, 104, 111, 118, 125]
    img_files = [data_dir + "obj{0:04d}_cc.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - d, 15s
    fnum = [89, 96, 103, 110, 117, 124]
    img_files = [data_dir + "obj{0:04d}_cd.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - e, 15s
    fnum = [88, 95, 102, 109, 116, 123]
    img_files = [data_dir + "obj{0:04d}_ce.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')
    
    # TT, 15s
    fnum = [92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_r_15s.fits', flat_frame = flat_dir + 'flat.fits')

    return
    
    
def find_stars_pleiades():
    data_dir = root_dir + 'reduce/pleiades/'

    # Open Loop
    fnum = [31, 38, 45, 52, 59, 66, 94, 101, 108, 115, 129]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - a
    fnum = [30, 37, 44, 51, 58, 65, 72, 93, 100, 107, 114, 128]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - b
    fnum = [28, 35, 42, 49, 56, 63, 70, 91, 98, 105, 112, 119]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - c
    fnum = [27, 34, 41, 48, 55, 62, 69, 90, 97, 104, 111, 118, 125]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - d
    fnum = [26, 33, 40, 47, 54, 61, 68, 89, 96, 103, 110, 117, 124]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - e
    fnum = [25, 32, 39, 46, 53, 60, 67, 88, 95, 102, 109, 116, 123]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    # TT
    fnum = [29, 36, 43, 50, 57, 64, 71, 92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

    util.mkdir(stats_dir)
    
    ##########
    # R-band
    ##########
    # Open Loop
    fnum = [31, 38, 45, 52, 59, 66]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop - a
    fnum = [30, 37, 44, 51, 58, 65, 72]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closeda.fits')

    # Closed Loop - b
    fnum = [28, 35, 42, 49, 56, 63, 70]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedb.fits')

    # Closed Loop - c
    fnum = [27, 34, 41, 48, 55, 62, 69]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedc.fits')

    # Closed Loop - d
    fnum = [26, 33, 40, 47, 54, 61, 68]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedd.fits')

    # Closed Loop - e
    fnum = [25, 32, 39, 46, 53, 60, 67]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closede.fits')

    # TT 
    fnum = [29, 36, 43, 50, 57, 64, 71]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedtt.fits')

    ##########
    # I-band
    ##########
    # Open Loop
    fnum = [94, 101, 108, 115, 129]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop - a
    fnum = [93, 100, 107, 114, 128]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closeda.fits')

    # Closed Loop - b
    fnum = [91, 98, 105, 112, 119]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedb.fits')

    # Closed Loop - c
    fnum = [90, 97, 104, 111, 118, 125]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedc.fits')

    # Closed Loop - d
    fnum = [89, 96, 103, 110, 117, 124]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedd.fits')

    # Closed Loop - e
    fnum = [88, 95, 102, 109, 116, 123]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closede.fits')

    # TT 
    fnum = [92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedtt.fits')
    
    return



def stack_pleiades_band():
    """
    Stack all the images for the Open loop, TTF closed loop, and GLAO closed loop (a)
    for both R and I band data. 
    """
    data_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

    # Open loop R-band
    open_img_nums = [31, 38, 45, 52, 59, 66]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in open_img_nums]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in open_img_nums]
    open_output_root = data_dir + 'pleiades_stack_open_r'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop - a
    closed_img_nums = [30, 37, 44, 51, 58, 65, 72]
    closed_images = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in closed_img_nums]
    closed_starlists = [data_dir + 'obj{0:04d}_ca_clean_stars.txt'.format(ii) for ii in closed_img_nums]
    closed_output_root = data_dir + 'pleiades_stack_closed_r'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # TT 
    tt_img_nums = [29, 36, 43, 50, 57, 64, 71]
    tt_images = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in tt_img_nums]
    tt_starlists = [data_dir + "obj{0:04d}_tt_clean_stars.txt".format(ii) for ii in tt_img_nums]
    tt_output_root = data_dir + 'pleiades_stack_tt_r'
    reduce_fli.shift_and_add(tt_images, tt_starlists, tt_output_root, method='mean')

    ##########
    # I-band
    ##########
    # Open Loop
    open_img_nums_i = [94, 101, 108, 115, 129]
    open_images_i = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in open_img_nums_i]
    open_starlists_i = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in open_img_nums_i]
    open_output_root_i = data_dir + 'pleiades_stack_open_i'
    reduce_fli.shift_and_add(open_images_i, open_starlists_i, open_output_root_i, method='mean')

    # Closed Loop - a
    closed_img_nums_i = [93, 100, 107, 114, 128]
    closed_images_i = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in closed_img_nums_i]
    closed_starlists_i = [data_dir + 'obj{0:04d}_ca_clean_stars.txt'.format(ii) for ii in closed_img_nums_i]
    closed_output_root_i = data_dir + 'pleiades_stack_closed_i'
    reduce_fli.shift_and_add(closed_images_i, closed_starlists_i, closed_output_root_i, method='mean')

    # TT 
    tt_img_nums_i = [92, 99, 106, 113, 120, 127]
    tt_images_i = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in tt_img_nums_i]
    tt_starlists_i = [data_dir + "obj{0:04d}_tt_clean_stars.txt".format(ii) for ii in tt_img_nums_i]
    tt_output_root_i = data_dir + 'pleiades_stack_tt_i'
    reduce_fli.shift_and_add(tt_images_i, tt_starlists_i, tt_output_root_i, method='mean')

    return
    
def analyze_stacks():
    data_dir = root_dir + 'reduce/pleiades/'
    
    img_files = [data_dir + 'pleiades_stack_open_r.fits',
                 data_dir + 'pleiades_stack_tt_r.fits',
                 data_dir + 'pleiades_stack_closed_r.fits',
                 data_dir + 'pleiades_stack_open_i.fits',
                 data_dir + 'pleiades_stack_tt_i.fits',
                 data_dir + 'pleiades_stack_closed_i.fits']
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats='stats_stacks.fits')

    return
    
