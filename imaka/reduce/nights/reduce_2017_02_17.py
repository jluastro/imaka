import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
import reduce_fli
import calib
import util
import pdb
import os

root_dir = '/Users/fatimaabdurrahman/Desktop/20170217/FLI/'


def make_flat():
    flat_raw_dir = root_dir + 'twilights/'
    dark_raw_dir = root_dir + 'darks/'
    flat_out_dir = root_dir + "reduce/calib/"

    #util.mkdir(flat_out_dir)
    
    flat_num = np.arange(9, 21)
    flat_frames = ['{0:s}twi_{1:04d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:04d}.fits'.format(dark_raw_dir, ss) for ss in flat_num]
    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat.fits')

    return


def make_sky():
    data_dir = root_dir + 'Pleiades/'
    sky_dir = root_dir + 'reduce/sky/'

    #util.mkdir(out_dir)

    sky_num = [93, 94, 95, 96, 97, 154, 155, 156, 157, 158]
    sky_num += [242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252]
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky.fits')
    
    return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    #util.mkdir(out_dir)

    # Open Loop
    fnum = [22, 28, 34, 40, 46]
    fnum += [52, 58, 65, 72, 79, 86]
    fnum += [98, 105, 112, 119, 126, 133, 140, 147]
    fnum += [159, 166, 173, 180, 187, 195, 202, 209]
    fnum += [216, 223, 230, 237]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop - a
    fnum = [27, 33, 39, 45, 51]
    fnum += [57, 64, 71, 78, 85]
    fnum += [104, 111, 118, 125, 132, 139, 146]
    fnum += [153, 165, 172, 179, 194, 201, 208]
    fnum += [215, 222, 229, 236]
    img_files = [data_dir + 'obj{0:04d}_ca.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - b
    fnum = [25, 31, 37, 43, 49]
    fnum += [55, 62, 69, 76, 83, 90]
    fnum += [102, 109, 116, 123, 130, 137, 144, 151]
    fnum += [163, 170, 177, 184, 191, 199, 206, 213]
    fnum += [220, 227, 234, 241]
    img_files = [data_dir + "obj{0:04d}_cb.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - c
    fnum = [24, 30, 36, 42, 48]
    fnum += [54, 61, 68, 75, 82, 89]
    fnum += [101, 108, 115, 122, 129, 136, 143, 150]
    fnum += [162, 169, 176, 183, 190, 198, 205, 212]
    fnum += [219, 226, 233, 240]
    img_files = [data_dir + "obj{0:04d}_cc.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - d
    fnum = [23, 29, 35, 41, 47]
    fnum += [53, 60, 67, 74, 81, 88]
    fnum += [100, 107, 114, 121, 128, 135, 142, 149]
    fnum += [161, 168, 175, 182, 189, 197, 204, 211]
    fnum += [218, 225, 232, 239]
    img_files = [data_dir + "obj{0:04d}_cd.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - f
    fnum = [26, 32, 38, 44, 50]
    fnum += [56, 63, 70, 77, 84, 91, 92]
    fnum += [103, 110, 117, 124, 131, 138, 145, 152]
    fnum += [164, 171, 178, 185, 192, 200, 207, 214]
    fnum += [221, 228, 235]
    img_files = [data_dir + "obj{0:04d}_cf.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')
    

    return
    
    
def find_stars_pleiades():
    data_dir = root_dir + 'reduce/pleiades/'

    # Open Loop
    fnum = [22, 28, 34, 40, 46]
    fnum += [52, 58, 65, 72, 79, 86]
    fnum += [98, 105, 112, 119, 126, 133, 140, 147]
    fnum += [159, 166, 173, 180, 187, 195, 202, 209]
    fnum += [216, 223, 230, 237]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - a
    fnum = [27, 33, 39, 45, 51]
    fnum += [57, 64, 71, 78, 85]
    fnum += [104, 111, 118, 125, 132, 139, 146]
    fnum += [153, 165, 172, 179, 194, 201, 208]
    fnum += [215, 222, 229, 236]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - b
    fnum = [25, 31, 37, 43, 49]
    fnum += [55, 62, 69, 76, 83, 90]
    fnum += [102, 109, 116, 123, 130, 137, 144, 151]
    fnum += [163, 170, 177, 184, 191, 199, 206, 213]
    fnum += [220, 227, 234, 241]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - c
    fnum = [24, 30, 36, 42, 48]
    fnum += [54, 61, 68, 75, 82, 89]
    fnum += [101, 108, 115, 122, 129, 136, 143, 150]
    fnum += [162, 169, 176, 183, 190, 198, 205, 212]
    fnum += [219, 226, 233, 240]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - d
    fnum = [23, 29, 35, 41, 47]
    fnum += [53, 60, 67, 74, 81, 88]
    fnum += [100, 107, 114, 121, 128, 135, 142, 149]
    fnum += [161, 168, 175, 182, 189, 197, 204, 211]
    fnum += [218, 225, 232, 239]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - f
    fnum = [26, 32, 38, 44, 50]
    fnum += [56, 63, 70, 77, 84, 91, 92]
    fnum += [103, 110, 117, 124, 131, 138, 145, 152]
    fnum += [164, 171, 178, 185, 192, 200, 207, 214]
    fnum += [221, 228, 235]
    img_files = [data_dir + "obj{0:04d}_cf_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
    fnum = [22, 28, 34, 40, 46]
    fnum += [52, 58, 65, 72, 79, 86]
    fnum += [98, 105, 112, 119, 126, 133, 140, 147]
    fnum += [159, 166, 173, 180, 187, 195, 202, 209]
    fnum += [216, 223, 230, 237]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop - a
    fnum = [27, 33, 39, 45, 51]
    fnum += [57, 64, 71, 78, 85]
    fnum += [104, 111, 118, 125, 132, 139, 146]
    fnum += [153, 165, 172, 179, 194, 201, 208]
    fnum += [215, 222, 229, 236]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closeda.fits')

    # Closed Loop - b
    fnum = [25, 31, 37, 43, 49]
    fnum += [55, 62, 69, 76, 83, 90]
    fnum += [102, 109, 116, 123, 130, 137, 144, 151]
    fnum += [163, 170, 177, 184, 191, 199, 206, 213]
    fnum += [220, 227, 234, 241]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedb.fits')

    # Closed Loop - c
    fnum = [24, 30, 36, 42, 48]
    fnum += [54, 61, 68, 75, 82, 89]
    fnum += [101, 108, 115, 122, 129, 136, 143, 150]
    fnum += [162, 169, 176, 183, 190, 198, 205, 212]
    fnum += [219, 226, 233, 240]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedc.fits')

    # Closed Loop - d
    fnum = [23, 29, 35, 41, 47]
    fnum += [53, 60, 67, 74, 81, 88]
    fnum += [100, 107, 114, 121, 128, 135, 142, 149]
    fnum += [161, 168, 175, 182, 189, 197, 204, 211]
    fnum += [218, 225, 232, 239]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedd.fits')

    # Closed Loop - f
    fnum = [26, 32, 38, 44, 50]
    fnum += [56, 63, 70, 77, 84, 91, 92]
    fnum += [103, 110, 117, 124, 131, 138, 145, 152]
    fnum += [164, 171, 178, 185, 192, 200, 207, 214]
    fnum += [221, 228, 235]
    img_files = [data_dir + "obj{0:04d}_cf_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedf.fits')

    
    return