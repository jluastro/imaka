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

root_dir = '/Users/fatimaabdurrahman/Desktop/20170218/FLI/'


def make_sky():
    data_dir = root_dir + 'Pleiades/'
    sky_dir = root_dir + 'reduce/sky/'

    #util.mkdir(out_dir)

    # 15 sec integration time
    sky_num = np.arange(73, 83)
    sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_dir+'pleiades_sky_15s.fits')
    
    return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    #util.mkdir(out_dir)

#     # Open Loop 15 s
#     fnum = [31, 38, 45, 52, 59, 66, 94, 101, 108, 115, 129]
#     img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame=flat_dir + 'flat.fits')
    
#     # Closed Loop - a, 15s
#     fnum = [30, 37, 44, 51, 58, 65, 72, 93, 100, 107, 114, 128]
#     img_files = [data_dir + 'obj{0:04d}_ca.fits'.format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - b, 15s
    fnum = [28, 35, 42, 49, 56, 63, 70, 91, 98, 105, 112, 119]
    img_files = [data_dir + "obj{0:04d}_cb.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - c, 15s
    fnum = [27, 34, 41, 48, 55, 62, 69, 90, 97, 104, 111, 118, 125]
    img_files = [data_dir + "obj{0:04d}_cc.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - d, 15s
    fnum = [26, 33, 40, 47, 54, 61, 68, 89, 96, 103, 110, 117, 124]
    img_files = [data_dir + "obj{0:04d}_cd.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed Loop - e, 15s
    fnum = [25, 32, 39, 46, 53, 60, 67, 88, 95, 102, 109, 116, 123]
    img_files = [data_dir + "obj{0:04d}_ce.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')
    
    # TT, 15s
    fnum = [29, 36, 43, 50, 57, 64, 71, 92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_15s.fits', flat_frame = flat_dir + 'flat.fits')
    

    return
    
    
def find_stars_pleiades():
    data_dir = root_dir + 'reduce/pleiades/'

#     # Open Loop
#     fnum = [31, 38, 45, 52, 59, 66, 94, 101, 108, 115, 129]
#     img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

#     # Closed Loop - a
#     fnum = [30, 37, 44, 51, 58, 65, 72, 93, 100, 107, 114, 128]
#     img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

#     # Closed Loop - b
#     fnum = [28, 35, 42, 49, 56, 63, 70, 91, 98, 105, 112, 119]
#     img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

#     # Closed Loop - c
#     fnum = [27, 34, 41, 48, 55, 62, 69, 90, 97, 104, 111, 118, 125]
#     img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

#     # Closed Loop - d
#     fnum = [26, 33, 40, 47, 54, 61, 68, 89, 96, 103, 110, 117, 124]
#     img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

#     # Closed Loop - e
#     fnum = [25, 32, 39, 46, 53, 60, 67, 88, 95, 102, 109, 116, 123]
#     img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    # TT
    fnum = [29, 36, 43, 50, 57, 64, 71, 92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
#     fnum = [31, 38, 45, 52, 59, 66, 94, 101, 108, 115, 129]
#     img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

#     # Closed Loop - a
#     fnum = [30, 37, 44, 51, 58, 65, 72, 93, 100, 107, 114, 128]
#     img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closeda.fits')

#     # Closed Loop - b
#     fnum = [28, 35, 42, 49, 56, 63, 70, 91, 98, 105, 112, 119]
#     img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedb.fits')

#     # Closed Loop - c
#     fnum = [27, 34, 41, 48, 55, 62, 69, 90, 97, 104, 111, 118, 125]
#     img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedc.fits')

#     # Closed Loop - d
#     fnum = [26, 33, 40, 47, 54, 61, 68, 89, 96, 103, 110, 117, 124]
#     img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
#     reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedd.fits')

    # Closed Loop - e
    fnum = [25, 32, 39, 46, 53, 60, 67, 88, 95, 102, 109, 116, 123]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closede.fits')

    # TT 
    fnum = [29, 36, 43, 50, 57, 64, 71, 92, 99, 106, 113, 120, 127]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedtt.fits')

    return